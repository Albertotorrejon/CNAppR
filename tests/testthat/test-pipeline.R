test_that("end-to-end pipeline runs without errors", {
  skip_on_cran()

  test_data <- data.frame(
    ID = rep(c("sample_1", "sample_2"), each = 10),
    chr = rep(1:10, 2),
    loc.start = rep(seq(1000000, 100000000, length.out = 10), 2),
    loc.end   = rep(seq(2000000, 101000000, length.out = 10), 2),
    seg.mean  = c(rnorm(10, 0.5, 0.05), rnorm(10, -0.3, 0.05)),
    BAF       = rep(0.5, 20),
    purity    = rep(0.8, 20),
    stringsAsFactors = FALSE
  )

  validation <- validate_cna_data(test_data)
  expect_true(validation$valid)

  var_prep <- prepare_clinical_variables(test_data)
  expect_type(var_prep, "list")
  expect_true("sample_names" %in% names(var_prep))

  reseg <- resegment_sample(test_data, sample_id = "sample_1")
  expect_s3_class(reseg, "data.frame")
  expect_true("classified" %in% colnames(reseg))

  reseg_list <- list(
    sample_1 = resegment_sample(test_data, sample_id = "sample_1"),
    sample_2 = resegment_sample(test_data, sample_id = "sample_2")
  )
  scores <- calculate_cna_scores(reseg_list)
  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 2)
})

test_that("validate_cna_data() catches missing required columns", {
  bad_data <- data.frame(ID = 1:5, chr = 1:5)
  validation <- validate_cna_data(bad_data)
  expect_false(validation$valid)
  expect_gt(length(validation$errors), 0)
})

test_that(
  "get_cytobands_data() returns correct structure for level3 and level4", {
  l3 <- get_cytobands_data("level3")
  expect_s3_class(l3, "data.frame")
  expect_true(all(c("chr", "label", "length") %in% colnames(l3)))
  expect_equal(nrow(l3), 48) # 2 arms x 24 chromosomes

  l4 <- get_cytobands_data("level4")
  expect_s3_class(l4, "data.frame")
  expect_equal(nrow(l4), 24) # one row per chromosome
  expect_true("length" %in% colnames(l4))
})


# ─── .load_targets_bed() ─────────────────────────────────────────────────────

test_that(
  ".load_targets_bed() parses a valid BED file and converts coordinates", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  tmp <- tempfile(fileext = ".bed")
  writeLines(
    c("chr1\t1000\t2000", "chr1\t5000\t6000", "chr2\t100000\t200000"),
    tmp
  )
  on.exit(unlink(tmp))

  chr_lens <- c(chr1 = 249250621L, chr2 = 243199373L)
  result   <- CNAppR:::.load_targets_bed(tmp, chr_lens)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 3)
  # BED 0-based → GRanges 1-based: start must be +1
  expect_equal(GenomicRanges::start(result)[1], 1001L)
  expect_equal(GenomicRanges::end(result)[1],   2000L)
})

test_that(".load_targets_bed() adds chr prefix when BAM uses it", {
  skip_if_not_installed("GenomicRanges")

  tmp <- tempfile(fileext = ".bed")
  writeLines(c("1\t1000\t2000", "2\t5000\t6000"), tmp)
  on.exit(unlink(tmp))

  chr_lens <- c(chr1 = 249250621L, chr2 = 243199373L)
  result   <- CNAppR:::.load_targets_bed(tmp, chr_lens)

  expect_equal(length(result), 2)
  expect_true(all(grepl("^chr", as.character(GenomicRanges::seqnames(result)))))
})

test_that(".load_targets_bed() strips chr prefix when BAM does not use it", {
  skip_if_not_installed("GenomicRanges")

  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t1000\t2000", "chr2\t5000\t6000"), tmp)
  on.exit(unlink(tmp))

  chr_lens <- c("1" = 249250621L, "2" = 243199373L)
  result   <- CNAppR:::.load_targets_bed(tmp, chr_lens)

  expect_equal(length(result), 2)
  seqnames_chr <- as.character(GenomicRanges::seqnames(result))
  expect_false(any(grepl("^chr", seqnames_chr)))
})

test_that(".load_targets_bed() skips comment, track and browser lines", {
  skip_if_not_installed("GenomicRanges")

  tmp <- tempfile(fileext = ".bed")
  writeLines(
    c("# comment", "track name=targets", "browser position chr1:1-1000",
      "chr1\t1000\t2000"),
    tmp
  )
  on.exit(unlink(tmp))

  chr_lens <- c(chr1 = 249250621L)
  result   <- CNAppR:::.load_targets_bed(tmp, chr_lens)

  expect_equal(length(result), 1)
})

test_that(".load_targets_bed() silently drops chromosomes not in the BAM", {
  skip_if_not_installed("GenomicRanges")

  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t1000\t2000", "chrUn_random\t500\t1000"), tmp)
  on.exit(unlink(tmp))

  chr_lens <- c(chr1 = 249250621L)
  result   <- CNAppR:::.load_targets_bed(tmp, chr_lens)

  expect_equal(length(result), 1)
  expect_equal(as.character(GenomicRanges::seqnames(result)), "chr1")
})

test_that(".load_targets_bed() errors on missing file", {
  expect_error(
    CNAppR:::.load_targets_bed("nonexistent_file.bed", c(chr1 = 249250621L)),
    "not found"
  )
})

test_that(".load_targets_bed() errors when no targets match BAM chromosomes", {
  tmp <- tempfile(fileext = ".bed")
  writeLines("chrUn\t1000\t2000", tmp)
  on.exit(unlink(tmp))

  expect_error(
    CNAppR:::.load_targets_bed(tmp, c(chr1 = 249250621L, chr2 = 243199373L)),
    "No BED targets"
  )
})


# ─── WES vs WGS log2_ratio behaviour ────────────────────────────────────────

test_that("read_bam_qc() WES: bins with 0 reads receive NA, not a large negative log2", {
  # Integration test: requires a real BAM + BED.
  # Expected behaviour: in WES mode (targets_bed != NULL), bins with
  # reads == 0 must have log2_ratio == NA instead of log2(0.5/median),
  # so CBS ignores them rather than treating them as deep deletions.
  skip("Requires a real indexed BAM and a capture kit BED file")
})

test_that("read_bam_qc() min_reads filters bins with marginal coverage", {
  # Integration test: requires a real BAM.
  # Expected behaviour: bins with reads < min_reads must be set to NA.
  skip("Requires a real indexed BAM file")
})
