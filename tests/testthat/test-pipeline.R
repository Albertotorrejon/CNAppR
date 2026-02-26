context("Full Pipeline")

test_that("end-to-end pipeline runs without errors", {
  skip_on_cran()

  # Create minimal test data
  test_data <- data.frame(
    ID = rep(c("sample_1", "sample_2"), each = 10),
    chr = rep(1:10, 2),
    loc.start = rep(seq(1000, 100000, length.out = 10), 2),
    loc.end = rep(seq(2000, 101000, length.out = 10), 2),
    seg.mean = c(rnorm(10, 0.5), rnorm(10, -0.3)),
    BAF = rep(0.5, 20),
    purity = rep(0.8, 20),
    stringsAsFactors = FALSE
  )

  # Validate data
  validation <- validate_cna_data(test_data)
  expect_true(validation$valid)

  # Prepare variables
  var_prep <- prepare_clinical_variables(test_data)
  expect_is(var_prep, "list")
  expect_true("sample_names" %in% names(var_prep))

  # Resegment sample
  reseg <- resegment_sample(test_data, sample_id = "sample_1")
  expect_is(reseg, "data.frame")
  expect_true("classified" %in% colnames(reseg))

  # Calculate scores
  reseg_list <- list(
    resegment_sample(test_data, sample_id = "sample_1"),
    resegment_sample(test_data, sample_id = "sample_2")
  )
  names(reseg_list) <- c("sample_1", "sample_2")

  scores <- calculate_cna_scores(reseg_list)
  expect_is(scores, "data.frame")
  expect_equal(nrow(scores), 2)
})

test_that("data validation catches errors", {
  bad_data <- data.frame(
    ID = 1:5,
    chr = 1:5
    # Missing required columns
  )

  validation <- validate_cna_data(bad_data)
  expect_false(validation$valid)
  expect_true(length(validation$errors) > 0)
})
