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

test_that("get_cytobands_data() returns correct structure for level3 and level4", {
  l3 <- get_cytobands_data("level3")
  expect_s3_class(l3, "data.frame")
  expect_true(all(c("chr", "label", "length") %in% colnames(l3)))
  expect_equal(nrow(l3), 48) # 2 arms x 24 chromosomes

  l4 <- get_cytobands_data("level4")
  expect_s3_class(l4, "data.frame")
  expect_equal(nrow(l4), 24) # one row per chromosome
  expect_true("length" %in% colnames(l4))
})
