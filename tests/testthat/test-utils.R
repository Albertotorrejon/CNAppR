test_that("read_cna_file() reads a valid TSV correctly", {
  tmp <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(ID = "s1", chr = 1, loc.start = 1000, loc.end = 2000,
               seg.mean = 0.5, stringsAsFactors = FALSE),
    tmp, sep = "\t", row.names = FALSE, quote = FALSE
  )
  df <- read_cna_file(tmp)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_true("seg.mean" %in% colnames(df))
  unlink(tmp)
})

test_that("read_cna_file() errors on missing file", {
  expect_error(read_cna_file("nonexistent_file.tsv"))
})

test_that("validate_cna_data() passes with correct data", {
  df <- data.frame(ID = "s1", chr = 1, loc.start = 1000, loc.end = 2000,
                   seg.mean = 0.5, stringsAsFactors = FALSE)
  result <- validate_cna_data(df)
  expect_true(result$valid)
  expect_length(result$errors, 0)
})

test_that("validate_cna_data() catches missing required columns", {
  bad <- data.frame(ID = 1:5, chr = 1:5)
  result <- validate_cna_data(bad)
  expect_false(result$valid)
  expect_gt(length(result$errors), 0)
})

test_that("validate_cna_data() catches non-numeric seg.mean", {
  df <- data.frame(ID = "s1", chr = 1, loc.start = 1000, loc.end = 2000,
                   seg.mean = "high", stringsAsFactors = FALSE)
  result <- validate_cna_data(df)
  expect_false(result$valid)
})

test_that("get_cytobands_data() returns correct structure for level3", {
  l3 <- get_cytobands_data("level3")
  expect_s3_class(l3, "data.frame")
  expect_true(all(c("chr", "label", "start", "end", "length") %in% colnames(l3)))
  expect_equal(nrow(l3), 48)
})

test_that("get_cytobands_data() returns correct structure for level4", {
  l4 <- get_cytobands_data("level4")
  expect_s3_class(l4, "data.frame")
  expect_true(all(c("chr", "start", "end", "length") %in% colnames(l4)))
  expect_equal(nrow(l4), 24)
})

test_that("get_cytobands_data() errors on unknown level", {
  expect_error(get_cytobands_data("level99"))
})
