# Tests for CNAppR_segments S3 class, .estimate_coverage, and run_bam_pipeline

# ─── CNAppR_segments S3 class ────────────────────────────────────────────────

test_that("CNAppR_segments constructor and print work", {
  segs <- data.frame(ID = "S01", chr = 1L, loc.start = 1L,
                     loc.end = 1e6L, seg.mean = 0.1)
  obj <- CNAppR:::new_cnapp_segments(
    segments = segs, method = "CBS", coverage = 35.2
  )
  expect_s3_class(obj, "CNAppR_segments")
  expect_equal(obj$method, "CBS")
  expect_equal(obj$coverage, 35.2)
  expect_true(is.na(obj$tumor_fraction))
  expect_true(is.na(obj$ploidy))
  expect_output(print(obj), "CNAppR_segments")
  expect_output(print(obj), "CBS")
})

test_that("as.data.frame.CNAppR_segments returns segment table", {
  segs <- data.frame(ID = "S01", chr = 1L, loc.start = 1L,
                     loc.end = 1e6L, seg.mean = 0.1)
  obj <- CNAppR:::new_cnapp_segments(segs, "CBS", 30)
  df  <- as.data.frame(obj)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1L)
  expect_named(df, names(segs))
})

test_that("CNAppR_segments with ichorCNA stores TF and ploidy", {
  segs <- data.frame(ID = "LB01", chr = 1L, loc.start = 1L,
                     loc.end = 1e6L, seg.mean = 0.05)
  obj <- CNAppR:::new_cnapp_segments(
    segments = segs, method = "ichorCNA", coverage = 0.8,
    tumor_fraction = 0.12, ploidy = 1.95
  )
  expect_equal(obj$tumor_fraction, 0.12)
  expect_equal(obj$ploidy, 1.95)
  expect_output(print(obj), "Tumor fraction")
  expect_output(print(obj), "ichorCNA")
})

test_that("CBS result has NA tumor_fraction and ploidy", {
  segs <- data.frame(ID = "S01", chr = 1L, loc.start = 1L,
                     loc.end = 1e8L, seg.mean = 0)
  obj  <- CNAppR:::new_cnapp_segments(segs, "CBS", 45)
  expect_true(is.na(obj$tumor_fraction))
  expect_true(is.na(obj$ploidy))
})

test_that("ichorCNA result has non-NA tumor_fraction and ploidy", {
  segs <- data.frame(ID = "LB01", chr = 1L, loc.start = 1L,
                     loc.end = 1e8L, seg.mean = 0.03)
  obj  <- CNAppR:::new_cnapp_segments(segs, "ichorCNA", 1.2,
                                       tumor_fraction = 0.07, ploidy = 2.1)
  expect_false(is.na(obj$tumor_fraction))
  expect_false(is.na(obj$ploidy))
})


# ─── .estimate_coverage ──────────────────────────────────────────────────────

test_that(".estimate_coverage returns numeric >= 0", {
  qc <- data.frame(
    loc.start = seq(1, by = 500000L, length.out = 10),
    loc.end   = seq(500000L, by = 500000L, length.out = 10),
    reads     = c(150, 160, 145, 170, 155, 165, 140, 158, 162, 148)
  )
  cov <- CNAppR:::.estimate_coverage(qc, read_length = 150L)
  expect_true(is.numeric(cov))
  expect_true(cov > 0)
})

test_that(".estimate_coverage returns 0 for all-NA reads", {
  qc <- data.frame(
    loc.start = 1:5, loc.end = 500000L * 1:5, reads = rep(NA_real_, 5)
  )
  expect_equal(CNAppR:::.estimate_coverage(qc), 0)
})

test_that(".estimate_coverage returns 0 for all-zero reads", {
  qc <- data.frame(
    loc.start = 1:5, loc.end = 500000L * 1:5, reads = rep(0, 5)
  )
  expect_equal(CNAppR:::.estimate_coverage(qc), 0)
})


# ─── run_bam_pipeline routing ────────────────────────────────────────────────

test_that("run_bam_pipeline routes CBS when coverage >= threshold", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("DNAcopy")
  skip_if_not_installed("mockery")

  mockery::stub(CNAppR::run_bam_pipeline, "read_bam_qc", function(...) {
    n <- 50L
    data.frame(
      chr        = rep("1", n),
      loc.start  = seq(1, by = 1e6L, length.out = n),
      loc.end    = seq(1e6L, by = 1e6L, length.out = n),
      reads      = rep(200L, n),
      log2_ratio = rnorm(n, 0, 0.05),
      bin_type   = "wgs"
    )
  })
  mockery::stub(CNAppR::run_bam_pipeline, ".estimate_coverage",
                function(...) 35)
  mockery::stub(CNAppR::run_bam_pipeline, "run_cbs_segmentation",
                function(...) {
                  data.frame(ID = "S01", chr = 1L,
                             loc.start = 1L, loc.end = 1e8L, seg.mean = 0)
                })
  mockery::stub(CNAppR::run_bam_pipeline, "resegment_sample",
                function(data, ...) data)

  res <- CNAppR::run_bam_pipeline("fake.bam", sample_id = "S01",
                                   method = "auto",
                                   coverage_threshold = 10)
  expect_s3_class(res, "CNAppR_segments")
  expect_equal(res$method, "CBS")
  expect_true(is.na(res$tumor_fraction))
  expect_true(is.na(res$ploidy))
})

test_that("run_bam_pipeline routes ichorCNA when coverage < threshold", {
  skip_if_not_installed("HMMcopy")
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("mockery")

  mockery::stub(CNAppR::run_bam_pipeline, "read_bam_qc", function(...) {
    n <- 100L
    data.frame(
      chr        = rep(as.character(1:22), length.out = n),
      loc.start  = seq(1, by = 5e5L, length.out = n),
      loc.end    = seq(5e5L, by = 5e5L, length.out = n),
      reads      = rep(5L, n),
      log2_ratio = rnorm(n, 0, 0.1),
      bin_type   = "wgs"
    )
  })
  mockery::stub(CNAppR::run_bam_pipeline, ".estimate_coverage",
                function(...) 0.5)
  mockery::stub(CNAppR::run_bam_pipeline, ".run_ichorcna", function(...) {
    segs <- data.frame(ID = "LB01", chr = 1L,
                       loc.start = 1L, loc.end = 1e8L, seg.mean = 0.02)
    CNAppR:::new_cnapp_segments(segs, "ichorCNA", 0.5, 0.08, 2.0)
  })

  res <- CNAppR::run_bam_pipeline("fake_ctdna.bam", sample_id = "LB01",
                                   method = "auto",
                                   coverage_threshold = 10)
  expect_s3_class(res, "CNAppR_segments")
  expect_equal(res$method, "ichorCNA")
  expect_false(is.na(res$tumor_fraction))
  expect_false(is.na(res$ploidy))
})

test_that("run_bam_pipeline forces CBS when method='CBS' regardless of coverage", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("DNAcopy")
  skip_if_not_installed("mockery")

  mockery::stub(CNAppR::run_bam_pipeline, "read_bam_qc", function(...) {
    data.frame(chr = "1", loc.start = 1L, loc.end = 5e5L,
               reads = 3L, log2_ratio = 0, bin_type = "wgs")
  })
  mockery::stub(CNAppR::run_bam_pipeline, ".estimate_coverage",
                function(...) 0.1)
  mockery::stub(CNAppR::run_bam_pipeline, "run_cbs_segmentation",
                function(...) {
                  data.frame(ID = "S01", chr = 1L,
                             loc.start = 1L, loc.end = 1e8L, seg.mean = 0)
                })
  mockery::stub(CNAppR::run_bam_pipeline, "resegment_sample",
                function(data, ...) data)

  res <- CNAppR::run_bam_pipeline("fake.bam", sample_id = "S01",
                                   method = "CBS",
                                   coverage_threshold = 10)
  expect_equal(res$method, "CBS")
})

test_that("run_bam_pipeline forces ichorCNA when method='ichorCNA'", {
  skip_if_not_installed("HMMcopy")
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("mockery")

  mockery::stub(CNAppR::run_bam_pipeline, "read_bam_qc", function(...) {
    data.frame(chr = "1", loc.start = 1L, loc.end = 5e5L,
               reads = 500L, log2_ratio = 0.02, bin_type = "wgs")
  })
  mockery::stub(CNAppR::run_bam_pipeline, ".estimate_coverage",
                function(...) 50)
  mockery::stub(CNAppR::run_bam_pipeline, ".run_ichorcna", function(...) {
    segs <- data.frame(ID = "S01", chr = 1L,
                       loc.start = 1L, loc.end = 1e8L, seg.mean = 0)
    CNAppR:::new_cnapp_segments(segs, "ichorCNA", 50, 0.3, 2.1)
  })

  res <- CNAppR::run_bam_pipeline("fake.bam", sample_id = "S01",
                                   method = "ichorCNA",
                                   coverage_threshold = 10)
  expect_equal(res$method, "ichorCNA")
})

test_that("run_bam_pipeline routes CBS when coverage equals threshold", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("DNAcopy")
  skip_if_not_installed("mockery")

  mockery::stub(CNAppR::run_bam_pipeline, "read_bam_qc", function(...) {
    data.frame(chr = "1", loc.start = 1L, loc.end = 5e5L,
               reads = 50L, log2_ratio = 0, bin_type = "wgs")
  })
  mockery::stub(CNAppR::run_bam_pipeline, ".estimate_coverage",
                function(...) 10)
  mockery::stub(CNAppR::run_bam_pipeline, "run_cbs_segmentation",
                function(...) {
                  data.frame(ID = "S01", chr = 1L,
                             loc.start = 1L, loc.end = 1e8L, seg.mean = 0)
                })
  mockery::stub(CNAppR::run_bam_pipeline, "resegment_sample",
                function(data, ...) data)

  res <- CNAppR::run_bam_pipeline("fake.bam", sample_id = "S01",
                                   method = "auto",
                                   coverage_threshold = 10)
  expect_equal(res$method, "CBS")
})
