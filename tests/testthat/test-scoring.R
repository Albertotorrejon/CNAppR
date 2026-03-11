test_that("calculate_cna_scores() returns correct structure", {
  reseg1 <- data.frame(
    ID = "sample_1", chr = 1, loc.start = 1000, loc.end = 2000,
    seg.mean = 0.5, classified = "Focal", stringsAsFactors = FALSE
  )
  reseg2 <- data.frame(
    ID = "sample_2", chr = 1, loc.start = 1000, loc.end = 2000,
    seg.mean = -0.3, classified = "Arm", stringsAsFactors = FALSE
  )
  reseg_list <- list(sample_1 = reseg1, sample_2 = reseg2)
  scores <- calculate_cna_scores(reseg_list)

  expect_s3_class(scores, "data.frame")
  expect_equal(ncol(scores), 3)
  expect_equal(nrow(scores), 2)
  expect_true(all(c("FCS", "BCS", "GCS") %in% colnames(scores)))
  expect_equal(rownames(scores), c("sample_1", "sample_2"))
})

test_that("classify_cna() correctly classifies gains and losses", {
  args <- list(0.2, log2(3/2), 1.0, -0.2, log2(1/2), log2(0.6/2))

  expect_equal(do.call(classify_cna, c(list(0.25), args)), 1)  # Low gain
  expect_equal(do.call(classify_cna, c(list(0.8),  args)), 2)  # Medium gain
  expect_equal(do.call(classify_cna, c(list(1.2),  args)), 3)  # High gain
  expect_equal(do.call(classify_cna, c(list(-0.25), args)), 1) # Low loss
  expect_equal(do.call(classify_cna, c(list(0.05), args)), 0)  # No alteration
})

test_that("calculate_cna_scores() handles mismatched sample_names length", {
  reseg_list <- list(s1 = data.frame(classified = "Focal"))
  expect_error(calculate_cna_scores(reseg_list, sample_names = c("a", "b")))
})
