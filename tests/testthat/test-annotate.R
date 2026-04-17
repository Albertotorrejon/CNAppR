test_that("annotate_cna_genes() adds genes column to output", {
  skip_on_cran()

  seg <- data.frame(
    chr = 19, loc.start = 42930402, loc.end = 45033565,
    classified = "focal", type = "Gain",
    stringsAsFactors = FALSE
  )

  result <- annotate_cna_genes(seg, genome_build = "hg19", only_cna = TRUE)

  expect_s3_class(result, "data.frame")
  expect_true("genes" %in% colnames(result))
  expect_false(is.na(result$genes[1]))
})

test_that("annotate_cna_genes() returns NA genes for diploid segments when only_cna = TRUE", {
  skip_on_cran()

  seg <- data.frame(
    chr = 1, loc.start = 1000000, loc.end = 2000000,
    classified = "diploid", type = "diploid",
    stringsAsFactors = FALSE
  )

  result <- annotate_cna_genes(seg, genome_build = "hg19", only_cna = TRUE)

  expect_true("genes" %in% colnames(result))
  expect_true(is.na(result$genes[1]))
})

test_that("annotate_cna_genes() errors on invalid genome_build", {
  seg <- data.frame(
    chr = 1, loc.start = 1000000, loc.end = 2000000,
    classified = "focal", type = "Gain",
    stringsAsFactors = FALSE
  )
  expect_error(annotate_cna_genes(seg, genome_build = "hg99"))
})
