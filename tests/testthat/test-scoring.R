context("CNA Scoring")

test_that("calculate_cna_scores() returns correct structure", {
  # Create mock resegmented data
  reseg1 <- data.frame(
    ID = "sample_1",
    chr = 1,
    loc.start = 1000,
    loc.end = 2000,
    seg.mean = 0.5,
    classified = "Focal",
    stringsAsFactors = FALSE
  )

  reseg2 <- data.frame(
    ID = "sample_2",
    chr = 1,
    loc.start = 1000,
    loc.end = 2000,
    seg.mean = -0.3,
    classified = "Arm",
    stringsAsFactors = FALSE
  )

  reseg_list <- list(sample_1 = reseg1, sample_2 = reseg2)

  scores <- calculate_cna_scores(reseg_list)

  expect_is(scores, "data.frame")
  expect_equal(ncol(scores), 3) # FCS, BCS, GCS
  expect_equal(nrow(scores), 2)
  expect_true(all(c("FCS", "BCS", "GCS") %in% colnames(scores)))
  expect_equal(rownames(scores), c("sample_1", "sample_2"))
})

test_that("classify_cna() correctly classifies gains and losses", {
  # Test gain classification
  w_low <- classify_cna(0.25, 0.2, log2(3/2), 1.0, -0.2, log2(1/2), log2(0.6/2))
  expect_equal(w_low, 1) # Low gain

  w_med <- classify_cna(0.8, 0.2, log2(3/2), 1.0, -0.2, log2(1/2), log2(0.6/2))
  expect_equal(w_med, 2) # Medium gain

  w_high <- classify_cna(1.2, 0.2, log2(3/2), 1.0, -0.2, log2(1/2), log2(0.6/2))
  expect_equal(w_high, 3) # High gain

  # Test loss classification
  w_loss_low <- classify_cna(-0.25, 0.2, log2(3/2), 1.0, -0.2, log2(1/2), log2(0.6/2))
  expect_equal(w_loss_low, 1) # Low loss

  # Test no alteration
  w_none <- classify_cna(0.05, 0.2, log2(3/2), 1.0, -0.2, log2(1/2), log2(0.6/2))
  expect_equal(w_none, 0) # No alteration
})
