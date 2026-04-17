test_that("calculate_cna_scores() returns correct structure", {
  reseg1 <- data.frame(classified = c("focal", "diploid"), score = c(4, 0),
                       stringsAsFactors = FALSE)
  reseg2 <- data.frame(classified = c("arm", "diploid"), score = c(2, 0),
                       stringsAsFactors = FALSE)
  scores <- calculate_cna_scores(list(sample_1 = reseg1, sample_2 = reseg2))

  expect_s3_class(scores, "data.frame")
  expect_equal(ncol(scores), 3)
  expect_equal(nrow(scores), 2)
  expect_true(all(c("FCS", "BCS", "GCS") %in% colnames(scores)))
  expect_equal(rownames(scores), c("sample_1", "sample_2"))
})

test_that("calculate_cna_scores() computes FCS and BCS correctly", {
  reseg1 <- data.frame(classified = c("focal", "focal"), score = c(2, 3),
                       stringsAsFactors = FALSE)
  reseg2 <- data.frame(classified = c("arm", "chromosomal"), score = c(1, 2),
                       stringsAsFactors = FALSE)
  scores <- calculate_cna_scores(list(s1 = reseg1, s2 = reseg2))

  expect_equal(scores["s1", "FCS"], 5)
  expect_equal(scores["s1", "BCS"], 0)
  expect_equal(scores["s2", "FCS"], 0)
  expect_equal(scores["s2", "BCS"], 3)
})

test_that("classify_cna() correctly classifies gains and losses", {
  low_g <- 0.2; med_g <- 0.58; high_g <- 1.0
  low_l <- -0.2; med_l <- -1.0; high_l <- -1.74

  expect_equal(classify_cna( 0.05, low_g, med_g, high_g, low_l, med_l, high_l), 0)
  expect_equal(classify_cna( 0.25, low_g, med_g, high_g, low_l, med_l, high_l), 1)
  expect_equal(classify_cna( 0.80, low_g, med_g, high_g, low_l, med_l, high_l), 2)
  expect_equal(classify_cna( 1.20, low_g, med_g, high_g, low_l, med_l, high_l), 3)
  expect_equal(classify_cna(-0.25, low_g, med_g, high_g, low_l, med_l, high_l), 1)
  expect_equal(classify_cna(-1.20, low_g, med_g, high_g, low_l, med_l, high_l), 2)
  expect_equal(classify_cna(-2.00, low_g, med_g, high_g, low_l, med_l, high_l), 3)
})

test_that("calculate_cna_scores() handles mismatched sample_names length", {
  reseg_list <- list(s1 = data.frame(classified = "focal", score = 1,
                                     stringsAsFactors = FALSE))
  expect_error(calculate_cna_scores(reseg_list, sample_names = c("a", "b")))
})
