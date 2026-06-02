# ichorcna.R — IchorCNA integration for low-coverage / liquid biopsy samples
#
# Adapted from IchorCNA v0.2.0 (Adalsteinsson et al., Nature Communications 2017)
# Original authors: Gavin Ha, Justin Rhoades — Dana-Farber / Broad Institute
# Source: https://github.com/broadinstitute/ichorCNA
# Modifications: functions ported to operate on plain data frames instead of
# RangedData objects; renamed with .ichor_ prefix to avoid namespace conflicts.
#
# Hard dependency: HMMcopy (Bioconductor) — required for the C-compiled Viterbi
# and forward-backward algorithms (.Call("viterbi", PACKAGE="HMMcopy"), etc.).
# Install with: BiocManager::install("HMMcopy")


# ─── pure math helpers (unchanged from IchorCNA) ─────────────────────────────

.ichor_normalize <- function(A) {
  vn <- function(x) x / (sum(x) + (sum(x) == 0))
  if (length(dim(A)) < 2) vn(A) else t(apply(A, 1, vn))
}

.ichor_dirichletpdflog <- function(x, k) {
  lgamma(sum(k, na.rm = TRUE)) - sum(lgamma(k), na.rm = TRUE) +
    sum((k - 1) * log(x), na.rm = TRUE)
}

.ichor_gammapdflog <- function(x, a, b) {
  a * log(b) - lgamma(a) + (a - 1) * log(x) + (-b * x)
}

.ichor_betapdflog <- function(x, a, b) {
  -lbeta(a, b) + (a - 1) * log(x) + (b - 1) * log(1 - x)
}

.ichor_tdistPDF <- function(x, mu, lambda, nu) {
  S <- ncol(x)
  if (!is.null(S)) {
    p <- NULL
    for (s in seq_len(S)) {
      tpdf <- (gamma(nu / 2 + 0.5) / gamma(nu / 2)) *
        ((lambda[s] / (pi * nu)) ^ 0.5) *
        (1 + (lambda[s] * (x[, s] - mu[s]) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
      p <- cbind(p, tpdf)
    }
  } else {
    nu <- rep(nu, length(mu))
    p  <- (gamma(nu / 2 + 0.5) / gamma(nu / 2)) *
      ((lambda / (pi * nu)) ^ 0.5) *
      (1 + (lambda * (x - mu) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
  }
  p[is.na(p)] <- 1
  p
}

.ichor_get2and3ComponentMixture <- function(ct, ct.sc.status, n, sp, phi) {
  S  <- length(n)
  cn <- 2
  mu <- matrix(NA, nrow = nrow(ct), ncol = ncol(ct))
  for (s in seq_len(S)) {
    sc <- ct.sc.status[, s] == TRUE
    mu[sc, s]  <- log(((1 - n[s]) * (1 - sp[s]) * ct[sc, s] +
                         (1 - n[s]) * sp[s] * cn + n[s] * cn) /
                        ((1 - n[s]) * phi[s] + n[s] * cn))
    mu[!sc, s] <- log(((1 - n[s]) * ct[!sc, s] + n[s] * cn) /
                        ((1 - n[s]) * phi[s] + n[s] * cn))
  }
  mu
}

.ichor_stLikelihood <- function(n, sp, phi, lambda, params, D, rho) {
  KS  <- nrow(rho)
  lik <- 0
  lambdaKS <- as.matrix(expand.grid(as.data.frame(lambda)))
  mus      <- as.matrix(
    .ichor_get2and3ComponentMixture(params$jointCNstates,
                                     params$jointSCstatus, n, sp, phi)
  )
  for (ks in seq_len(KS)) {
    probs <- log(.ichor_tdistPDF(D, mus[ks, ], lambdaKS[ks, ], params$nu))
    lik   <- lik + as.numeric(rho[ks, ] %*% rowSums(as.data.frame(probs)))
  }
  lik
}

.ichor_priorProbs <- function(n, sp, phi, lambda, piG, A, params,
                               estimateNormal     = TRUE,
                               estimatePloidy     = TRUE,
                               estimatePrecision  = TRUE,
                               estimateTransition = TRUE,
                               estimateInitDist   = TRUE,
                               estimateSubclone   = TRUE) {
  S  <- params$numberSamples
  K  <- length(params$ct)
  KS <- K ^ S
  priorA <- 0
  if (estimateTransition)
    for (ks in seq_len(KS))
      priorA <- priorA + .ichor_dirichletpdflog(A[ks, ], params$dirPrior[ks, ])
  priorLambda <- 0
  if (estimatePrecision)
    for (s in seq_len(S))
      for (k in seq_len(K))
        priorLambda <- priorLambda +
          .ichor_gammapdflog(as.data.frame(lambda)[k, s],
                              params$alphaLambda[k], params$betaLambda[k, s])
  priorN   <- if (estimateNormal)  sum(.ichor_betapdflog(n,  params$alphaN, params$betaN))  else 0
  priorSP  <- if (estimateSubclone) sum(.ichor_betapdflog(sp, params$alphaSp, params$betaSp)) else 0
  priorPhi <- if (estimatePloidy)  sum(.ichor_gammapdflog(phi, params$alphaPhi, params$betaPhi)) else 0
  priorPi  <- if (estimateInitDist) .ichor_dirichletpdflog(piG, params$kappa) else 0
  list(prior       = priorA + as.numeric(priorLambda) + priorN + priorSP + priorPhi + priorPi,
       priorLambda = as.numeric(priorLambda))
}

.ichor_completeLikelihoodFun <- function(x, pType, n, sp, phi, lambda, piG, A,
                                          params, D, rho, Zcounts,
                                          estimateNormal     = TRUE,
                                          estimatePloidy     = TRUE,
                                          estimatePrecision  = TRUE,
                                          estimateTransition = FALSE,
                                          estimateInitDist   = TRUE,
                                          estimateSubclone   = TRUE) {
  KS     <- nrow(rho)
  S      <- params$numberSamples
  lambda <- matrix(lambda, ncol = S, byrow = FALSE)
  if (estimatePrecision  && any(pType == "lambda")) lambda <- matrix(x[pType == "lambda"], ncol = S, byrow = FALSE)
  if (estimateNormal     && any(pType == "n"))      n      <- x[pType == "n"]
  if (estimateSubclone   && any(pType == "sp"))     sp     <- x[pType == "sp"]
  if (estimatePloidy     && any(pType == "phi"))    phi    <- x[pType == "phi"]

  prior   <- .ichor_priorProbs(n, sp, phi, lambda, piG, A, params,
                                 estimateNormal, estimatePloidy, estimatePrecision,
                                 estimateTransition, estimateInitDist, estimateSubclone)
  likObs  <- .ichor_stLikelihood(n, sp, phi, lambda, params, D, rho)
  likInit <- rho[, 1] %*% log(piG)
  likTxn  <- 0
  for (ks in seq_len(KS))
    likTxn <- likTxn + Zcounts[ks, ] %*% log(A[ks, ])
  as.numeric(likObs + likInit + likTxn + prior$prior)
}

.ichor_estimateMixWeightsParamMap <- function(rho, kappa) {
  (rowSums(rho, na.rm = TRUE) + kappa - 1) /
    (sum(rowSums(rho, na.rm = TRUE)) + sum(kappa) - length(kappa))
}

.ichor_estimateParamsMap <- function(D, n_prev, sp_prev, phi_prev, lambda_prev,
                                      pi_prev, A_prev, params, rho, Zcounts,
                                      estimateNormal     = TRUE,
                                      estimatePloidy     = TRUE,
                                      estimatePrecision  = TRUE,
                                      estimateInitDist   = TRUE,
                                      estimateTransition = TRUE,
                                      estimateSubclone   = TRUE) {
  KS  <- nrow(rho)
  K   <- length(params$ct)
  S   <- length(n_prev)
  iN  <- c(1e-6, 1 - 1e-6)
  iP  <- c(.Machine$double.eps, 10)
  iL  <- c(1e-5, 1e4)

  n_hat      <- n_prev; sp_hat <- sp_prev; phi_hat <- phi_prev
  lambda_hat <- lambda_prev; pi_hat <- pi_prev; A_hat <- A_prev

  if (estimateTransition)
    for (k in seq_len(KS)) {
      A_hat[k, ] <- .ichor_normalize(Zcounts[k, ] + params$dirPrior[k, ])
    }
  if (estimateInitDist)
    pi_hat <- .ichor_estimateMixWeightsParamMap(rho, params$kappa)

  opt_args <- list(params = params, D = D, rho = rho, Zcounts = Zcounts,
                   estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                   estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                   estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                   method = "L-BFGS-B", control = list(trace = 0, fnscale = -1))

  if (estimateNormal) {
    est <- suppressWarnings(
      do.call(stats::optim, c(list(par = n_prev, fn = .ichor_completeLikelihoodFun,
                                   pType = rep("n", S), n = n_prev, sp = sp_prev,
                                   phi = phi_prev, lambda = lambda_prev, piG = pi_hat,
                                   A = A_hat, lower = iN[1], upper = iN[2]), opt_args))
    )
    n_hat <- est$par
  }
  if (estimateSubclone) {
    est <- suppressWarnings(
      do.call(stats::optim, c(list(par = sp_prev, fn = .ichor_completeLikelihoodFun,
                                   pType = rep("sp", S), n = n_hat, sp = sp_prev,
                                   phi = phi_prev, lambda = lambda_prev, piG = pi_hat,
                                   A = A_hat, lower = iN[1], upper = iN[2]), opt_args))
    )
    sp_hat <- est$par
  }
  if (estimatePloidy) {
    est <- suppressWarnings(
      do.call(stats::optim, c(list(par = phi_prev, fn = .ichor_completeLikelihoodFun,
                                   pType = rep("phi", length(phi_prev)), n = n_hat,
                                   sp = sp_hat, phi = phi_prev, lambda = lambda_prev,
                                   piG = pi_hat, A = A_hat, lower = iP[1], upper = iP[2]),
                              opt_args))
    )
    phi_hat <- est$par
  }
  estLambda <- suppressWarnings(
    do.call(stats::optim, c(list(par = c(lambda_prev), fn = .ichor_completeLikelihoodFun,
                                  pType = rep("lambda", K * S), n = n_hat, sp = sp_hat,
                                  phi = phi_hat, lambda = lambda_prev, piG = pi_hat,
                                  A = A_hat, lower = iL[1]), opt_args))
  )
  lambda_hat <- matrix(estLambda$par, ncol = S, byrow = FALSE)

  list(n = n_hat, sp = sp_hat, phi = phi_hat, lambda = lambda_hat,
       piG = pi_hat, A = A_hat, F = estLambda$value)
}


# ─── HMM parameters ──────────────────────────────────────────────────────────

.ichor_getTransitionMatrix <- function(K, e, strength) {
  A <- matrix((1 - e[1]) / (K - 1), K, K)
  diag(A) <- e[1]
  A  <- .ichor_normalize(A)
  list(A = A, dirPrior = A * strength[1])
}

.ichor_getDefaultParameters <- function(x, maxCN = 7, ct.sc = NULL,
                                         ploidy = 2, e = 0.9999999,
                                         e.sameState = 10, strength = 10000000,
                                         includeHOMD = FALSE, n_0 = 0.5) {
  ct    <- if (includeHOMD) 0:maxCN else 1:maxCN
  param <- list(
    strength = strength, e = e,
    ct = c(ct, ct.sc),
    ct.sc.status = c(rep(FALSE, length(ct)), rep(TRUE, length(ct.sc))),
    phi_0 = 2, alphaPhi = 4, betaPhi = 1.5,
    n_0 = n_0, alphaN = 2, betaN = 2,
    sp_0 = 0.5, alphaSp = 2, betaSp = 2,
    lambda = as.matrix(rep(100, length(ct) + length(ct.sc))),
    nu = 2.1,
    kappa = rep(75, length(ct)),
    alphaLambda = 5
  )
  K <- length(param$ct)
  S <- if (!is.null(dim(x))) ncol(x) else 1L
  param$numberSamples <- S
  betaLambdaVal <- if (S > 1)
    (apply(x, 2, stats::sd, na.rm = TRUE) / sqrt(K)) ^ 2
  else
    (stats::sd(x, na.rm = TRUE) / sqrt(K)) ^ 2
  param$betaLambda <- matrix(betaLambdaVal, ncol = S, nrow = K, byrow = TRUE)
  param$alphaLambda <- rep(param$alphaLambda, K)

  logR.var <- 1 / betaLambdaVal
  if (S > 1) {
    param$lambda <- matrix(logR.var, nrow = K, ncol = S, byrow = TRUE,
                            dimnames = list(c(), if (!is.null(colnames(x))) colnames(x) else NULL))
  } else {
    param$lambda <- matrix(logR.var, K)
    param$lambda[param$ct %in% c(1, 3)] <- logR.var
    param$lambda[param$ct >= 4]          <- logR.var / 5
    param$lambda[param$ct == max(param$ct)] <- logR.var / 15
    if (any(param$ct.sc.status)) param$lambda[param$ct.sc.status] <- logR.var / 10
  }
  param$jointCNstates <- expand.grid(rep(list(param$ct), S))
  param$jointSCstatus <- expand.grid(rep(list(param$ct.sc.status), S))
  colnames(param$jointCNstates) <- paste0("Sample.", seq_len(S))
  colnames(param$jointSCstatus) <- paste0("Sample.", seq_len(S))

  txn <- .ichor_getTransitionMatrix(K ^ S, e, strength)
  cnStateDiff <- apply(param$jointCNstates, 1,
                       function(r) abs(max(r) - min(r)))
  if (e.sameState > 0 && S > 1) {
    txn$A[, cnStateDiff == 0] <- txn$A[, cnStateDiff == 0] * e.sameState * K
    txn$A[, cnStateDiff >= 3] <- txn$A[, cnStateDiff >= 3] / e.sameState / K
  }
  diag(txn$A) <- e
  txn$A <- .ichor_normalize(txn$A)
  param$A        <- txn$A
  param$dirPrior <- txn$A * strength[1]
  if (any(param$ct.sc.status)) {
    param$A[, param$ct.sc.status]        <- param$A[, param$ct.sc.status] / 10
    param$A        <- .ichor_normalize(param$A)
    param$dirPrior[, param$ct.sc.status] <- param$dirPrior[, param$ct.sc.status] / 10
  }
  if (includeHOMD) {
    param$A[1, 2:K] <- param$A[1, 2:K] * 1e-5
    param$A[2:K, 1] <- param$A[2:K, 1] * 1e-5
    param$A[1, 1]   <- param$A[1, 1]   * 1e-5
    param$A        <- .ichor_normalize(param$A)
    param$dirPrior <- param$A * param$strength
  }
  param$kappa <- rep(75, K ^ S)
  param$kappa[cnStateDiff == 0] <- param$kappa[cnStateDiff == 0] + 125
  param$kappa[cnStateDiff >= 3] <- param$kappa[cnStateDiff >= 3] - 50
  param$kappa[rowSums(param$jointCNstates == 2) == S] <- 800
  param
}


# ─── EM ──────────────────────────────────────────────────────────────────────

.ichor_runEM <- function(copy, chr, chrTrain, param, maxiter, verbose = FALSE,
                          estimateNormal     = TRUE, estimatePloidy     = TRUE,
                          estimatePrecision  = TRUE, estimateTransition = TRUE,
                          estimateInitDist   = TRUE, estimateSubclone   = TRUE,
                          likChangeConvergence = 1e-3) {
  S   <- param$numberSamples
  K   <- length(param$ct)
  KS  <- K ^ S
  N   <- nrow(copy)
  rho <- matrix(0, KS, N)
  py  <- matrix(0, KS, N)
  mus     <- array(0, dim = c(KS, S, maxiter))
  lambdas <- array(0, dim = c(K,  S, maxiter))
  phi  <- matrix(NA, length(param$phi_0), maxiter)
  n    <- matrix(NA, S, maxiter)
  sp   <- matrix(NA, S, maxiter)
  piG  <- matrix(0,  KS, maxiter)
  A    <- .ichor_normalize(param$A)
  Zcounts  <- matrix(0, KS, KS)
  loglik   <- rep(0, maxiter)
  converged <- FALSE

  chrs  <- levels(chr)
  chrsI <- lapply(chrs, function(ch) which(chr == ch))

  i <- 1L
  piG[, i]     <- .ichor_normalize(param$kappa)
  n[, i]       <- param$n_0
  sp[, i]      <- param$sp_0
  phi[, i]     <- param$phi_0
  lambdas[, , i] <- param$lambda
  lambdasKS    <- as.matrix(expand.grid(as.data.frame(lambdas[, , i])))
  mus[, , i]   <- as.matrix(
    .ichor_get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus,
                                     n[, i], sp[, i], phi[, i])
  )
  for (ks in seq_len(KS)) {
    probs    <- .ichor_tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
    py[ks, ] <- apply(probs, 1, prod)
  }
  loglik[i] <- -Inf

  while (!converged && i < maxiter) {
    i <- i + 1L
    if (verbose) message("runEM iter ", i - 1, ": E-step")
    Zcounts[] <- 0
    for (j in seq_along(chrsI)) {
      I <- intersect(chrsI[[j]], which(chrTrain))
      if (length(I) > 0) {
        out <- .Call("forward_backward", piG[, i - 1], A, py[, I],
                     PACKAGE = "HMMcopy")
        rho[, I] <- out$rho
        loglik[i] <- loglik[i] + out$loglik
        Zcounts   <- Zcounts + t(colSums(aperm(out$xi, c(3, 2, 1))))
      }
    }
    if (verbose) message("runEM iter ", i - 1, ": M-step")
    out <- .ichor_estimateParamsMap(
      copy[chrTrain, , drop = FALSE],
      n[, i - 1], sp[, i - 1], phi[, i - 1], lambdas[, , i - 1],
      piG[, i - 1], A, param, rho[, chrTrain, drop = FALSE], Zcounts,
      estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
      estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
      estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone
    )
    n[, i]         <- out$n;   sp[, i]  <- out$sp;  phi[, i] <- out$phi
    lambdas[, , i] <- out$lambda;  piG[, i] <- out$piG;  A <- out$A

    lambdasKS  <- as.matrix(expand.grid(as.data.frame(lambdas[, , i])))
    mus[, , i] <- as.matrix(
      .ichor_get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus,
                                       n[, i], sp[, i], phi[, i])
    )
    for (ks in seq_len(KS)) {
      probs    <- .ichor_tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
      py[ks, ] <- apply(probs, 1, prod)
    }
    prior     <- .ichor_priorProbs(n[, i], sp[, i], phi[, i], lambdas[, , i],
                                    piG[, i], A, param,
                                    estimateNormal, estimatePloidy, estimatePrecision,
                                    estimateTransition, estimateInitDist, estimateSubclone)
    loglik[i] <- loglik[i] + prior$prior
    if (verbose) message("runEM iter ", i - 1, " loglik: ", round(loglik[i], 4))
    if ((abs(loglik[i] - loglik[i - 1]) / abs(loglik[i])) < likChangeConvergence)
      converged <- TRUE
    if (loglik[i] < loglik[i - 1]) { i <- i - 1L; converged <- TRUE }
  }
  if (converged)
    for (j in seq_along(chrsI)) {
      I   <- chrsI[[j]]
      out <- .Call("forward_backward", piG[, i], A, py[, I], PACKAGE = "HMMcopy")
      rho[, I] <- out$rho
    }

  list(n = n[, 1:i, drop = FALSE], sp = sp[, 1:i, drop = FALSE],
       phi = phi[, 1:i, drop = FALSE], mus = mus[, , 1:i, drop = FALSE],
       lambdas = lambdas[, , 1:i, drop = FALSE], pi = piG[, 1:i, drop = FALSE],
       A = A, loglik = loglik[1:i], rho = rho, param = param, py = py, iter = i)
}


# ─── Viterbi and segmentation ────────────────────────────────────────────────

.ichor_runViterbi <- function(convergedParams, chr) {
  chrs  <- levels(chr)
  chrsI <- lapply(chrs, function(ch) which(chr == ch))
  N     <- ncol(convergedParams$py)
  Z     <- rep(0L, N)
  piG   <- convergedParams$pi[, convergedParams$iter]
  A     <- convergedParams$A
  segs  <- vector("list", length(chrs))
  for (c in seq_along(chrsI)) {
    I      <- chrsI[[c]]
    out    <- .Call("viterbi", log(piG), log(A), log(convergedParams$py[, I]),
                    PACKAGE = "HMMcopy")
    Z[I]   <- out$path
    segs[[c]] <- out$seg
  }
  list(segs = segs, states = Z)
}

# Adapted from IchorCNA segmentData(): operates on plain data frames instead
# of RangedData objects. Expects each element of dataGR to have columns:
# space (factor), start, end, copy.
.ichor_segmentData <- function(dataGR, validInd, states, convergedParams) {
  includeHOMD <- any(convergedParams$param$ct == 0)
  ev_names    <- if (includeHOMD)
    c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP", paste0("HLAMP", 2:25))
  else
    c("HETD","NEUT","GAIN","AMP","HLAMP", paste0("HLAMP", 2:25))

  states    <- states[validInd]
  S         <- length(dataGR)
  colNames  <- c("space", "start", "end", "copy")
  segList   <- list()

  for (i in seq_len(S)) {
    id     <- names(dataGR)[i]
    dataIn <- dataGR[[i]][validInd, ]

    # RLE per chromosome (chromosome order follows factor levels)
    chr_levels <- levels(dataIn$space)
    chr_levels <- chr_levels[chr_levels %in% unique(as.character(dataIn$space))]
    rleResults <- t(sapply(chr_levels, function(ch) {
      ind <- as.character(dataIn$space) == ch
      rle(states[ind])
    }))

    rleLengths <- unlist(rleResults[, "lengths"])
    rleValues  <- unlist(rleResults[, "values"])
    sampleDF   <- as.data.frame(dataIn)[, colNames, drop = FALSE]
    numSegs    <- length(rleLengths)

    segs <- as.data.frame(
      matrix(NA, ncol = 7, nrow = numSegs,
             dimnames = list(NULL, c("chr","start","end","state","event",
                                     "median","copy.number")))
    )
    prevInd <- 0L
    for (j in seq_len(numSegs)) {
      s_start  <- prevInd + 1L
      s_end    <- prevInd + rleLengths[j]
      segDF    <- sampleDF[s_start:s_end, , drop = FALSE]
      prevInd  <- s_end
      numR     <- nrow(segDF)
      segs[j, "chr"]   <- as.character(segDF[1, "space"])
      segs[j, "start"] <- segDF[1, "start"]
      segs[j, "state"] <- rleValues[j]
      cn <- convergedParams$param$jointCNstates[rleValues[j], i]
      segs[j, "copy.number"] <- cn
      if (as.character(segDF[1, "space"]) == as.character(segDF[numR, "space"])) {
        segs[j, "end"]    <- segDF[numR, "end"]
        segs[j, "median"] <- round(stats::median(segDF$copy, na.rm = TRUE), 6)
        segs[j, "event"]  <- if (includeHOMD) ev_names[cn + 1] else ev_names[cn]
      }
    }
    segList[[id]] <- segs
  }
  segList
}

# Adapted from IchorCNA HMMsegment(): replaces space()/start()/end() calls on
# RangedData with direct data frame column access.
.ichor_HMMsegment <- function(x, validInd = NULL, dataType = "copy",
                                param      = NULL,
                                chrTrain   = as.character(1:22),
                                maxiter    = 50L,
                                estimateNormal     = TRUE,
                                estimatePloidy     = TRUE,
                                estimatePrecision  = TRUE,
                                estimateSubclone   = TRUE,
                                estimateTransition = TRUE,
                                estimateInitDist   = TRUE,
                                logTransform       = FALSE,
                                verbose            = FALSE) {
  chr     <- x[[1]]$space   # factor column — replaces space(x[[1]])
  dataMat <- as.matrix(as.data.frame(lapply(x, "[[", dataType)))

  if (logTransform) {
    dataMat <- apply(dataMat, 2,
                     function(col) log(col / stats::median(col, na.rm = TRUE)))
  } else {
    dataMat <- log(2 ^ dataMat)   # log2 → natural log
  }
  for (i in seq_along(x)) x[[i]]$copy <- dataMat[, i]

  chrInd <- chr %in% chrTrain
  if (!is.null(validInd)) chrInd <- chrInd & validInd

  if (is.null(param))
    param <- .ichor_getDefaultParameters(dataMat[chrInd, , drop = FALSE])

  convergedParams <- .ichor_runEM(
    dataMat, chr, chrInd, param, maxiter, verbose,
    estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
    estimateSubclone = estimateSubclone, estimatePrecision = estimatePrecision,
    estimateTransition = estimateTransition, estimateInitDist = estimateInitDist
  )

  viterbiResults <- .ichor_runViterbi(convergedParams, chr)

  segs <- .ichor_segmentData(x, validInd, viterbiResults$states, convergedParams)

  ev_names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP", paste0("HLAMP", 2:25))

  S <- length(x)
  cnaList <- list()
  for (s in seq_len(S)) {
    id          <- names(x)[s]
    copyNumber  <- convergedParams$param$jointCNstates[viterbiResults$states, s]
    sc_status   <- convergedParams$param$jointSCstatus[viterbiResults$states, s]
    cnaList[[id]] <- data.frame(
      sample          = as.character(id),
      chr             = as.character(x[[s]]$space),   # replaces space(x[[s]])
      start           = x[[s]]$start,                  # replaces start(x[[s]])
      end             = x[[s]]$end,                    # replaces end(x[[s]])
      copy.number     = copyNumber,
      event           = ev_names[pmin(copyNumber + 1, length(ev_names))],
      logR            = round(log2(exp(dataMat[, s])), 4),
      subclone.status = as.numeric(sc_status),
      stringsAsFactors = FALSE
    )
    segs[[s]]$median          <- log2(exp(segs[[s]]$median))
    segs[[s]]$subclone.status <-
      convergedParams$param$jointSCstatus[as.integer(segs[[s]]$state), s]
  }
  convergedParams$segs <- segs
  list(cna = cnaList, results = convergedParams, viterbiResults = viterbiResults)
}


# ─── GC correction (adapted from IchorCNA correctReadCounts) ─────────────────

# Adapted from IchorCNA utils.R correctReadCounts(). Operates on a data frame
# with columns: space, start, end, reads, gc, (optionally) map.
# Returns the data frame with added columns: cor.gc, cor.map, copy, valid, ideal.
.ichor_correctReadCounts <- function(df, chrNormalize = as.character(1:22),
                                      mappability = 0.9, samplesize = 50000L,
                                      verbose = FALSE) {
  if (is.null(df$reads) || is.null(df$gc))
    stop("Data frame must have 'reads' and 'gc' columns for GC correction.")

  chrInd   <- as.character(df$space) %in% as.character(chrNormalize)
  df$valid <- TRUE
  df$valid[df$reads <= 0 | is.na(df$reads) | df$gc < 0 | is.na(df$gc)] <- FALSE
  df$ideal <- TRUE

  routlier <- 0.01; doutlier <- 0.001
  range_r  <- stats::quantile(df$reads[df$valid & chrInd],
                               probs = c(0, 1 - routlier), na.rm = TRUE)
  domain_g <- stats::quantile(df$gc[df$valid & chrInd],
                               probs = c(doutlier, 1 - doutlier), na.rm = TRUE)
  has_map  <- !is.null(df$map) && any(!is.na(df$map))

  if (has_map) {
    df$ideal[!df$valid | df$map < mappability |
               df$reads <= range_r[1] | df$reads > range_r[2] |
               df$gc < domain_g[1]   | df$gc > domain_g[2]] <- FALSE
  } else {
    df$ideal[!df$valid | df$reads <= range_r[1] | df$reads > range_r[2] |
               df$gc < domain_g[1]   | df$gc > domain_g[2]] <- FALSE
  }

  if (verbose) message("IchorCNA GC correction: LOESS fitting...")
  set    <- which(df$ideal & chrInd)
  select <- set[sample.int(length(set), min(length(set), as.integer(samplesize)))]
  rough  <- stats::loess(df$reads[select] ~ df$gc[select], span = 0.03)
  gc_seq <- seq(0, 1, by = 0.001)
  final  <- stats::loess(stats::predict(rough, gc_seq) ~ gc_seq, span = 0.3)
  df$cor.gc <- df$reads / stats::predict(final, df$gc)

  if (has_map) {
    if (verbose) message("IchorCNA: mappability correction...")
    coutlier <- 0.01
    range_c  <- stats::quantile(df$cor.gc[df$valid & chrInd],
                                 probs = c(0, 1 - coutlier), na.rm = TRUE)
    set2     <- which(df$cor.gc < range_c[2] & chrInd)
    sel2     <- set2[sample.int(length(set2),
                                 min(length(set2), as.integer(samplesize)))]
    map_fn   <- stats::approxfun(stats::lowess(df$map[sel2], df$cor.gc[sel2]))
    df$cor.map <- df$cor.gc / map_fn(df$map)
  } else {
    df$cor.map <- df$cor.gc
  }

  df$copy              <- df$cor.map
  df$copy[df$copy <= 0 | is.na(df$copy)] <- NA_real_
  df$copy              <- log2(df$copy)
  df
}


# ─── CNAppR ↔ IchorCNA adapter ───────────────────────────────────────────────

# Parse a fixedStep WIG file into a data.frame with columns:
# chr (UCSC, e.g. "chr1"), start (integer), value (numeric).
.read_wig_fixedStep <- function(path) {
  lines     <- readLines(path, warn = FALSE)
  chr_name  <- NA_character_
  pos_start <- NA_integer_
  pos_step  <- NA_integer_
  n         <- length(lines)
  chr_vec   <- character(n)
  start_vec <- integer(n)
  val_vec   <- numeric(n)
  k <- 0L
  cur_pos   <- NA_integer_
  cur_step  <- NA_integer_
  cur_chr   <- NA_character_
  for (ln in lines) {
    ln <- trimws(ln)
    if (nchar(ln) == 0L || startsWith(ln, "#") || startsWith(ln, "track")) next
    if (startsWith(ln, "fixedStep")) {
      cur_chr   <- sub(".*chrom=(\\S+).*",  "\\1", ln)
      cur_pos   <- as.integer(sub(".*start=(\\d+).*",  "\\1", ln))
      cur_step  <- as.integer(sub(".*step=(\\d+).*",   "\\1", ln))
    } else {
      k <- k + 1L
      chr_vec[k]   <- cur_chr
      start_vec[k] <- cur_pos
      val_vec[k]   <- as.numeric(ln)
      cur_pos      <- cur_pos + cur_step
    }
  }
  data.frame(chr   = chr_vec[seq_len(k)],
             start = start_vec[seq_len(k)],
             value = val_vec[seq_len(k)],
             stringsAsFactors = FALSE)
}

# Locate a WIG file: first try CNAppR's own extdata, then the ichorCNA package.
.find_wig <- function(filename) {
  p <- system.file("extdata", filename, package = "CNAppR")
  if (nzchar(p) && file.exists(p)) return(p)
  p <- system.file("extdata", filename, package = "ichorCNA")
  if (nzchar(p) && file.exists(p)) return(p)
  NULL
}

# Converts read_bam_qc() output to the named list of data frames expected by
# .ichor_HMMsegment().  Uses IchorCNA's WIG-derived GC + mappability values so
# that bin masking (gc=-1, map<0.9) is identical to runIchorCNA.R, then runs
# IchorCNA's two-step LOESS correction on raw read counts.
.bam_bins_to_hmm_input <- function(qc_data, sample_id,
                                    genome_build = "hg38",
                                    bin_size     = 1000000L) {
  # Normalise chr names to NCBI style ("1"..."22","X","Y")
  chr_str <- gsub("^chr", "", as.character(qc_data$chr))
  chr_str[chr_str == "23"] <- "X"
  chr_str[chr_str == "24"] <- "Y"

  # Genomic sort: numeric autosomes first, then X, Y
  chr_numeric <- suppressWarnings(as.integer(chr_str))
  auto_levels  <- as.character(sort(unique(chr_numeric[!is.na(chr_numeric)])))
  sex_levels   <- sort(unique(chr_str[is.na(chr_numeric)]))
  all_levels   <- c(auto_levels, sex_levels)
  all_levels   <- all_levels[all_levels %in% unique(chr_str)]

  # Use raw read counts (pre-CNAppR LOESS) for IchorCNA's own correction
  raw_reads <- if (!is.null(qc_data$reads_raw)) as.numeric(qc_data$reads_raw)
               else as.numeric(qc_data$reads)

  df <- data.frame(
    space  = factor(chr_str, levels = all_levels),
    start  = as.integer(qc_data$loc.start),
    end    = as.integer(qc_data$loc.end),
    reads  = raw_reads,
    stringsAsFactors = FALSE
  )

  # ── Load IchorCNA WIG files for GC and mappability ─────────────────────────
  bin_kb   <- as.integer(round(bin_size / 1000L))
  gc_file  <- sprintf("gc_%s_%dkb.wig",  genome_build, bin_kb)
  map_file <- sprintf("map_%s_%dkb.wig", genome_build, bin_kb)
  gc_path  <- .find_wig(gc_file)
  map_path <- .find_wig(map_file)

  # UCSC chr keys for WIG matching (WIG uses "chr1", NCBI space uses "1")
  ucsc_key <- paste0("chr", as.character(df$space))

  if (!is.null(gc_path)) {
    message("  Loading IchorCNA GC WIG: ", basename(gc_path))
    gc_wig  <- .read_wig_fixedStep(gc_path)
    wig_key <- paste(gc_wig$chr, gc_wig$start)
    df_key  <- paste(ucsc_key,   df$start)
    m       <- match(df_key, wig_key)
    df$gc   <- gc_wig$value[m]   # NA where bin not in WIG
    n_masked  <- sum(!is.na(df$gc) & df$gc < 0)
    n_missing <- sum(is.na(df$gc))
    message("    GC: ", n_masked, " bins masked (gc=-1), ",
            n_missing, " bins without WIG match")
  } else {
    message("  WIG not found: ", gc_file, " — falling back to bins RDS gc values")
    df$gc <- if (!is.null(qc_data$gc)) as.numeric(qc_data$gc) else NULL
  }

  if (!is.null(map_path)) {
    message("  Loading IchorCNA map WIG: ", basename(map_path))
    map_wig <- .read_wig_fixedStep(map_path)
    wig_key <- paste(map_wig$chr, map_wig$start)
    df_key  <- paste(ucsc_key,    df$start)
    m       <- match(df_key, wig_key)
    df$map  <- map_wig$value[m]
    n_lowmap <- sum(!is.na(df$map) & df$map < 0.9)
    message("    Map: ", n_lowmap, " bins with mappability < 0.9")
  }

  # ── Apply IchorCNA's GC + mappability correction ────────────────────────────
  has_gc <- !is.null(df$gc) && any(!is.na(df$gc) & df$gc > 0)
  if (has_gc) {
    df <- .ichor_correctReadCounts(df, verbose = FALSE)
  } else {
    # No GC: fall back to log2_ratio already computed by read_bam_qc()
    df$copy  <- as.numeric(qc_data$log2_ratio)
    df$valid <- !is.na(df$copy) & !is.na(df$reads) & df$reads > 0
    df$ideal <- df$valid
  }

  # Sort by chromosome then position for consistent RLE in segmentData
  chr_order <- as.integer(factor(df$space, levels = levels(df$space)))
  df <- df[order(chr_order, df$start), ]
  df$space <- droplevels(df$space)

  stats::setNames(list(df), sample_id)
}


# ─── main IchorCNA runner ────────────────────────────────────────────────────

# Internal runner called by run_bam_pipeline().
# Returns a CNAppR_segments object.
.run_ichorcna <- function(qc_data, sample_id, coverage, genome_build, bin_size,
                           maxiter = 50L, estimateNormal = TRUE,
                           estimatePloidy = TRUE,
                           normal_init = c(0.2, 0.4, 0.5, 0.6, 0.8), ...) {
  .check_pkg("HMMcopy")
  requireNamespace("HMMcopy", quietly = TRUE)

  # Resolve bin_size: if NULL, infer from median bin width in qc_data
  if (is.null(bin_size) || length(bin_size) == 0L) {
    bin_size <- as.integer(
      round(stats::median(qc_data$loc.end - qc_data$loc.start, na.rm = TRUE))
    )
    message("  .run_ichorcna: bin_size inferred from qc_data median bin width: ",
            bin_size, " bp")
  }
  bin_size <- as.integer(bin_size)

  hmm_input <- .bam_bins_to_hmm_input(qc_data, sample_id,
                                        genome_build = genome_build,
                                        bin_size     = bin_size)

  # Mirror runIchorCNA.R's valid flag after filterByMappabilityScore:
  #   valid = reads>0 & gc>=0  (basic), with low-map bins (map<0.9) treated as absent.
  # IchorCNA physically removes low-map bins before passing to HMMsegment, so
  # 'valid' in its output already implies map>=0.9. We replicate this by
  # explicitly requiring map>=0.9 when a map WIG was loaded.
  df_hmm  <- hmm_input[[1]]
  has_map <- !is.null(df_hmm$map) && any(!is.na(df_hmm$map))
  if (has_map) {
    validInd <- df_hmm$valid &
                (is.na(df_hmm$map) | df_hmm$map >= 0.9) &
                !is.na(df_hmm$copy)
  } else {
    validInd <- df_hmm$valid & !is.na(df_hmm$copy)
  }

  n_validInd <- sum(validInd, na.rm = TRUE)
  message("  Valid bins for HMM training: ", n_validInd,
          " (valid & map>=0.9 & copy not NA)")
  if (n_validInd < 100L)
    stop("Too few valid bins (", n_validInd, ") for IchorCNA after map filter. ",
         "Check that gc_correction = TRUE and the BAM has sufficient coverage.")

  # Mirror runIchorCNA.R exactly:
  # - getDefaultParameters() receives log2 data (not natural log)
  # - lambda is recomputed after getDefaultParameters() with IchorCNA's exact formula
  # - alphaLambda = 3 (lambdaScaleHyperParam default in runIchorCNA.R)
  # - solution selected by highest loglik, but only among solutions with
  #   subclonal CNA fraction <= maxFracCNASubclone (0.7) and
  #   subclonal genome fraction <= maxFracGenomeSubclone (0.5)
  normal_inits <- normal_init
  message("  Running HMM-EM with ", length(normal_inits),
          " initialization(s) n_0=", paste(normal_inits, collapse=","),
          " (max ", maxiter, " iter each)...")

  # log2 data for getDefaultParameters (matches runIchorCNA.R line 249-250)
  logR_mat     <- as.matrix(as.data.frame(lapply(hmm_input, "[[", "copy")))
  chrInd_train <- (hmm_input[[1]]$space %in% as.character(1:22)) & validInd
  logR_train   <- logR_mat[chrInd_train, , drop = FALSE]

  all_res  <- list()
  all_ll   <- numeric(0)
  all_sc   <- numeric(0)   # fraction of subclonal CNA bins
  all_gsc  <- numeric(0)   # fraction of subclonal genome bins

  for (n0 in normal_inits) {
    # Build params with log2 training data (runIchorCNA.R line 250-251)
    param <- .ichor_getDefaultParameters(logR_train, maxCN = 7, n_0 = n0)

    # Recompute lambda exactly as runIchorCNA.R line 260:
    #   apply(logR, 2, sd) on all bins remaining after map filter (incl. chrX/Y).
    #   IchorCNA's logR is all post-mapFilter bins; na.rm skips gc=-1 (copy=NA) bins.
    #   We use validInd (all chrs, not just chrTrain) to match those remaining bins.
    logR_var <- 1 / ((stats::sd(logR_mat[validInd, , drop = FALSE], na.rm = TRUE) /
                        sqrt(length(param$ct))) ^ 2)
    param$lambda <- rep(logR_var, length(param$ct))
    param$lambda[param$ct %in% c(1, 3)] <- logR_var
    param$lambda[param$ct >= 4]          <- logR_var / 5
    param$lambda[param$ct == max(param$ct)] <- logR_var / 15
    if (any(param$ct.sc.status))
      param$lambda[param$ct.sc.status] <- logR_var / 10
    param$alphaLambda <- rep(3, length(param$ct))  # lambdaScaleHyperParam = 3

    res <- tryCatch(
      .ichor_HMMsegment(
        hmm_input,
        validInd       = validInd,
        param          = param,
        maxiter        = as.integer(maxiter),
        estimateNormal = estimateNormal,
        estimatePloidy = estimatePloidy,
        verbose        = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      iter <- res$results$iter
      ll   <- res$results$loglik[iter]
      # subclonal fractions (runIchorCNA.R lines 330-334)
      cna       <- res$cna[[sample_id]]
      sc_bins   <- sum(cna$subclone.status, na.rm = TRUE)
      tot_bins  <- nrow(cna)
      alt_bins  <- sum(cna$copy.number != 2L, na.rm = TRUE)
      frac_gsc  <- sc_bins / max(tot_bins, 1L)
      frac_csc  <- if (alt_bins > 0) sc_bins / alt_bins else 0
      message("    n_0 = ", n0,
              "  TF = ",    round(1 - res$results$n[1, iter], 4),
              "  loglik = ", round(ll, 2),
              "  sc_frac = ", round(frac_csc, 3))
      all_res  <- c(all_res,  list(res))
      all_ll   <- c(all_ll,  ll)
      all_sc   <- c(all_sc,  frac_csc)
      all_gsc  <- c(all_gsc, frac_gsc)
    }
  }
  if (length(all_res) == 0)
    stop("IchorCNA HMM failed for all initializations.")

  # Select solution: prefer those passing subclonal thresholds (runIchorCNA.R 352-362)
  pass <- which(all_sc <= 0.7 & all_gsc <= 0.5)
  idx  <- if (length(pass) > 0) pass[which.max(all_ll[pass])] else which.max(all_ll)
  hmm_res     <- all_res[[idx]]
  best_loglik <- all_ll[[idx]]
  message("  Selected: TF = ",
          round(1 - hmm_res$results$n[1, hmm_res$results$iter], 4),
          "  loglik = ", round(best_loglik, 2))

  # Extract converged tumor fraction and ploidy
  iter <- hmm_res$results$iter
  normal_frac    <- hmm_res$results$n[1, iter]
  tumor_fraction <- 1 - normal_frac
  ploidy         <- hmm_res$results$phi[1, iter]

  message("  Converged: tumor fraction = ", round(tumor_fraction, 4),
          ", ploidy = ", round(ploidy, 4))

  # Convert IchorCNA segments to CNAppR format
  segs_ichor <- hmm_res$results$segs[[sample_id]]
  if (is.null(segs_ichor) || nrow(segs_ichor) == 0)
    stop("IchorCNA produced no segments for sample '", sample_id, "'.")

  seg_cnapp <- data.frame(
    ID        = sample_id,
    chr       = .norm_chr_to_int(segs_ichor$chr),
    loc.start = as.integer(segs_ichor$start),
    loc.end   = as.integer(segs_ichor$end),
    seg.mean  = as.numeric(segs_ichor$median),
    stringsAsFactors = FALSE
  )
  seg_cnapp <- seg_cnapp[!is.na(seg_cnapp$chr), ]

  message("segment_cna [3/3]: Re-segmentation and classification...")
  reseg <- resegment_sample(
    seg_cnapp,
    sample_id   = sample_id,
    sample_type = "liquid_biopsy",
    purity      = tumor_fraction,
    ...
  )

  # Apply TF/ploidy correction to bin log2 ratios so they align with segment
  # means in plot_genome_wide_cna() — identical to IchorCNA's display correction:
  # corrected = raw_log2R - log2(((1-TF)*2 + TF*ploidy) / ploidy)
  qc_plot <- qc_data
  correction <- log2(((1 - tumor_fraction) * 2 + tumor_fraction * ploidy) / ploidy)
  qc_plot$log2_ratio <- qc_data$log2_ratio - correction

  new_cnapp_segments(
    segments       = reseg,
    method         = "ichorCNA",
    coverage       = coverage,
    tumor_fraction = tumor_fraction,
    ploidy         = ploidy,
    metadata       = list(sample_id    = sample_id,
                          genome_build = genome_build,
                          bin_size     = bin_size,
                          qc_data      = qc_plot)
  )
}
