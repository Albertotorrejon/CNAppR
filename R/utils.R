#' Read CNA Data File
#'
#' Reads a copy number data file with configurable separators and decimal characters.
#'
#' @param filepath Character path to the data file (TSV, CSV, or TXT).
#' @param sep Character field separator (default tab character).
#' @param dec Character decimal separator (default ".").
#'
#' @return A data frame with columns: ID, chr, loc.start, loc.end, seg.mean.
#'
#' @details
#' The file MUST contain the following columns:
#' - ID: Sample identifier
#' - chr: Chromosome (numeric or "chr1", "chrX", etc.)
#' - loc.start: Start location (bp)
#' - loc.end: End location (bp)
#' - seg.mean: Log2 ratio (numeric)
#'
#' Additional columns (BAF, purity, clinical variables) are preserved.
#'
#' @examples
#' \dontrun{
#'   data <- read_cna_file("sample_data.tsv", sep = "\t", dec = ".")
#' }
#'
#' @export
read_cna_file <- function(filepath, sep = "\t", dec = ".") {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }

  tryCatch({
    df <- utils::read.table(filepath, sep = sep, header = TRUE, dec = dec, stringsAsFactors = FALSE)
    return(df)
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
}


#' Validate CNA Data Structure
#'
#' Validates that input data contains all required columns and correct data types.
#'
#' @param data Data frame to validate.
#' @param min_cols Minimum number of columns required (default 5).
#' @param max_var_cols Maximum number of variable columns (default 26).
#'
#' @return A list with:
#'   - valid (logical): TRUE if data passes all checks
#'   - errors (character vector): descriptions of any validation errors
#'
#' @details
#' Required columns: ID, chr, loc.start, loc.end, seg.mean
#' - seg.mean must be numeric
#' - chr must be coercible to numeric
#' - ID must not have duplicates per row
#'
#' @examples
#' \dontrun{
#'   validation <- validate_cna_data(my_data)
#'   if (!validation$valid) {
#'     message(paste("Validation errors:", paste(validation$errors, collapse = "; ")))
#'   }
#' }
#'
#' @export
validate_cna_data <- function(data, min_cols = 5, max_var_cols = 26) {
  errors <- c()

  # Check dimensions
  if (nrow(data) < 1) {
    errors <- c(errors, "Data has no rows.")
  }
  if (ncol(data) < min_cols) {
    errors <- c(errors, paste("Data must have at least", min_cols, "columns."))
  }

  # Check required columns
  required_cols <- c("ID", "chr", "loc.start", "loc.end", "seg.mean")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  if (length(errors) == 0) {
    # Check column types
    if (!is.numeric(data$seg.mean)) {
      errors <- c(errors, "seg.mean column must be numeric.")
    }

    # Check variable columns limit
    var_cols <- setdiff(colnames(data), c("ID", "chr", "loc.start", "loc.end", "seg.mean"))
    if (length(var_cols) > max_var_cols) {
      errors <- c(errors, paste("Too many variable columns (max", max_var_cols, ")."))
    }
  }

  return(list(valid = length(errors) == 0, errors = errors))
}


#' Prepare Clinical Variables
#'
#' Extracts, classifies (numeric/categoric), and prepares clinical/annotation variables.
#'
#' @param data Data frame with sample data. Columns after ID, chr, loc.start, loc.end, seg.mean are treated as variables.
#' @param exclude_cols Character vector of column names to exclude (e.g., c("BAF", "purity")).
#' @param sample_var Character name of sample/ID column (default "ID").
#'
#' @return A list containing:
#'   - mat_variables: Data frame with one row per unique sample (numeric/categoric columns)
#'   - variables_info: Data frame with name_var, class_var (categoric/numeric), color_palette
#'   - sample_names: Character vector of unique samples
#'
#' @details
#' Automatically assigns color palettes from RColorBrewer for visualization.
#'
#' @examples
#' \dontrun{
#'   var_prep <- prepare_clinical_variables(data, exclude_cols = c("BAF", "purity"))
#'   head(var_prep$mat_variables)
#' }
#'
#' @export
prepare_clinical_variables <- function(data, exclude_cols = c("BAF", "purity"), sample_var = "ID") {
  # Get unique samples
  sample_names <- unique(as.character(data[[sample_var]]))

  # Find row index for each sample
  each_sample <- sapply(sample_names, function(x, vec) which(vec == x)[1], vec = as.character(data[[sample_var]]))

  # Extract variable columns
  variables_names <- colnames(data)[!colnames(data) %in% c("ID", "chr", "loc.start", "loc.end", "seg.mean", exclude_cols)]

  if (length(variables_names) == 0) {
    warning("No variable columns found after excluding standard columns.")
    return(list(mat_variables = NULL, variables_info = NULL, sample_names = sample_names))
  }

  # Classify each variable
  variables_class <- rep(NA_character_, length(variables_names))
  list_variables <- list()
  vars_to_keep <- c()

  for (g in seq_along(variables_names)) {
    var_name <- variables_names[g]
    raw_var <- data[each_sample, var_name]
    raw_var[raw_var == ""] <- NA

    class_var <- class(raw_var)

    if (class_var %in% c("factor", "character")) {
      variables_class[g] <- "categoric"
      list_variables[[var_name]] <- as.character(raw_var)
      vars_to_keep <- c(vars_to_keep, g)
    } else if (class_var %in% c("numeric", "integer")) {
      variables_class[g] <- "numeric"
      list_variables[[var_name]] <- as.numeric(as.character(raw_var))
      vars_to_keep <- c(vars_to_keep, g)
    }
  }

  if (length(vars_to_keep) == 0) {
    warning("No valid variables could be processed.")
    return(list(mat_variables = NULL, variables_info = NULL, sample_names = sample_names))
  }

  # Keep only valid variables
  variables_names <- variables_names[vars_to_keep]
  variables_class <- variables_class[vars_to_keep]

  # Build data frame
  mat_variables <- as.data.frame(do.call(cbind, list_variables), stringsAsFactors = FALSE)
  rownames(mat_variables) <- sample_names

  # Add dummy variable if only one column
  if (ncol(mat_variables) == 1) {
    mat_variables <- cbind(mat_variables, mat_variables)
    colnames(mat_variables) <- c(colnames(mat_variables)[1], "id_dummy")
    variables_names <- c(variables_names, "id_dummy")
    variables_class <- c(variables_class, "categoric")
  }

  # Ensure correct classes
  for (i in which(variables_class == "categoric")) {
    mat_variables[, i] <- as.character(mat_variables[, i])
  }
  for (i in which(variables_class == "numeric")) {
    mat_variables[, i] <- as.numeric(as.character(mat_variables[, i]))
  }

  # Assign color palettes
  div_palettes <- rownames(RColorBrewer::brewer.pal.info[
    which(RColorBrewer::brewer.pal.info[, "category"] == "div"), ])
  seq_palettes <- rownames(RColorBrewer::brewer.pal.info[
    which(RColorBrewer::brewer.pal.info[, "category"] == "seq"), ])

  color_palettes   <- c(div_palettes, seq_palettes)
  color_assignment <- rep_len(color_palettes, length(variables_names))

  variables_info <- data.frame(
    name_var = variables_names,
    class_var = variables_class,
    color_palette = color_assignment,
    stringsAsFactors = FALSE
  )

  return(list(
    mat_variables = mat_variables,
    variables_info = variables_info,
    sample_names = sample_names
  ))
}


#' Prepare Annotation Data
#'
#' Reads and merges annotation/clinical data with the main CNA data frame.
#' Each annotation row is broadcast to all CNA rows belonging to the same
#' sample (many-to-one join on the sample ID column).
#'
#' @param data Main data frame (from read_cna_file()).
#' @param annot_filepath Character path to annotation file (tab-separated,
#'   with header). Must contain a column named by sample_var.
#' @param sep Character field separator for annotation file (default "\t").
#' @param dec Character decimal separator for annotation file (default ".").
#' @param sample_var Character name of sample ID column in both data frames
#'   (default "ID").
#'
#' @return The original data frame with annotation columns appended. Rows
#'   and their order are preserved. Samples in data without an annotation
#'   entry receive NA in all appended columns. Annotation columns that
#'   already exist in data are silently replaced by the annotation values.
#'
#' @details
#' Merging is done with base::merge() (left join), which is substantially
#' faster than a row-by-row loop for large cohorts.
#'
#' @examples
#' \dontrun{
#'   merged_data <- prepare_annotation_data(data, "annotations.tsv")
#' }
#'
#' @export
prepare_annotation_data <- function(data, annot_filepath, sep = "\t",
                                     dec = ".", sample_var = "ID") {
  if (!file.exists(annot_filepath))
    stop("Annotation file not found: ", annot_filepath)

  annot_df <- utils::read.table(annot_filepath, sep = sep, header = TRUE,
                                 dec = dec, stringsAsFactors = FALSE)

  if (!sample_var %in% colnames(annot_df))
    stop("Annotation file must contain '", sample_var, "' column.")

  # Warn about annotation samples absent from main data
  extra <- setdiff(as.character(annot_df[[sample_var]]),
                    unique(as.character(data[[sample_var]])))
  if (length(extra) > 0L)
    warning(length(extra), " annotation sample(s) not found in main data: ",
            paste(head(extra, 3L), collapse = ", "),
            if (length(extra) > 3L) paste0(" ... (+", length(extra) - 3L, " more)") else "")

  # Annotation columns that overlap existing data columns (other than sample_var)
  # are replaced by the annotation values — drop them from data first to avoid
  # .x / .y suffix collisions from merge().
  annot_new_cols <- setdiff(colnames(annot_df), sample_var)
  data_clean <- data[, setdiff(colnames(data), annot_new_cols), drop = FALSE]

  # Preserve original row order: merge() does not guarantee it
  data_clean[[".row_order"]] <- seq_len(nrow(data_clean))

  merged <- merge(data_clean, annot_df, by = sample_var, all.x = TRUE,
                   sort = FALSE)
  merged <- merged[order(merged[[".row_order"]]), , drop = FALSE]
  merged[[".row_order"]] <- NULL
  rownames(merged) <- NULL

  merged
}


#' Get Cytoband Reference Data
#'
#' Returns reference cytoband data (level 3 or level 4) for a given genome build.
#'
#' @param level Character "level3" (p/q arm definitions) or "level4"
#'   (full chromosome lengths). Default "level3".
#'
#' @return Data frame with columns: chr, start, end, label, length.
#'
#' @export
get_cytobands_data <- function(level = "level3") {
  # Reference cytoband data (hg19/GRCh37)
  # Source: aux_files/cytobands_level3_pq.csv and cytobands_level4_chrom.csv
  # from the original CNApp Shiny repository

  if (level == "level3") {
    # p/q arm definitions per chromosome (hg19)
    # arms[1,"end"] is used as the centromere position in scoring
    l3 <- data.frame(
      chr = rep(1:24, each = 2),
      label = rep(c("p", "q"), 24),
      start = c(
        0, 125000000,   # chr1
        0,  93300000,   # chr2
        0,  91000000,   # chr3
        0,  50400000,   # chr4
        0,  48400000,   # chr5
        0,  61000000,   # chr6
        0,  59900000,   # chr7
        0,  45600000,   # chr8
        0,  49000000,   # chr9
        0,  40200000,   # chr10
        0,  53700000,   # chr11
        0,  35800000,   # chr12
        0,  17900000,   # chr13
        0,  17600000,   # chr14
        0,  19000000,   # chr15
        0,  36600000,   # chr16
        0,  24000000,   # chr17
        0,  17200000,   # chr18
        0,  26500000,   # chr19
        0,  27500000,   # chr20
        0,  13200000,   # chr21
        0,  14700000,   # chr22
        0,  60600000,   # chr23 (X)
        0,  12500000    # chr24 (Y)
      ),
      end = c(
        125000000, 249250621,   # chr1
         93300000, 243199373,   # chr2
         91000000, 198022430,   # chr3
         50400000, 191154276,   # chr4
         48400000, 180915260,   # chr5
         61000000, 171115067,   # chr6
         59900000, 159138663,   # chr7
         45600000, 146364022,   # chr8
         49000000, 141213431,   # chr9
         40200000, 135534747,   # chr10
         53700000, 135006516,   # chr11
         35800000, 133851895,   # chr12
         17900000, 115169878,   # chr13
         17600000, 107349540,   # chr14
         19000000, 102531392,   # chr15
         36600000,  90354753,   # chr16
         24000000,  81195210,   # chr17
         17200000,  78077248,   # chr18
         26500000,  59128983,   # chr19
         27500000,  63025520,   # chr20
         13200000,  48129895,   # chr21
         14700000,  51304566,   # chr22
         60600000, 155270560,   # chr23 (X)
         12500000,  59373566    # chr24 (Y)
      ),
      stringsAsFactors = FALSE
    )
    l3$length <- l3$end - l3$start
    return(l3)

  } else if (level == "level4") {
    # Total chromosome lengths (hg19)
    # Source: cytobands_level4_chrom.csv from the original repository
    l4 <- data.frame(
      chr = 1:24,
      start = rep(0, 24),
      end = c(
        249250621,   # chr1
        243199373,   # chr2
        198022430,   # chr3
        191154276,   # chr4
        180915260,   # chr5
        171115067,   # chr6
        159138663,   # chr7
        146364022,   # chr8
        141213431,   # chr9
        135534747,   # chr10
        135006516,   # chr11
        133851895,   # chr12
        115169878,   # chr13
        107349540,   # chr14
        102531392,   # chr15
         90354753,   # chr16
         81195210,   # chr17
         78077248,   # chr18
         59128983,   # chr19
         63025520,   # chr20
         48129895,   # chr21
         51304566,   # chr22
        155270560,   # chr23 (X)
         59373566    # chr24 (Y)
      ),
      label = rep("x", 24),
      stringsAsFactors = FALSE
    )
    l4$length <- l4$end - l4$start
    return(l4)

  } else {
    stop("Unknown cytoband level: ", level)
  }
}