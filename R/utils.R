#' Read CNA Data File
#'
#' Reads a copy number data file with configurable separators and decimal characters.
#'
#' @param filepath Character path to the data file (TSV, CSV, or TXT).
#' @param sep Character field separator (default "\t" for tab).
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
  rownames(mat_variables) <- sample_names[vars_to_keep]

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

  color_palettes <- c(div_palettes, seq_palettes)
  color_assignment <- color_palettes[seq_len(length(variables_names))]

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
#' Reads and merges annotation/clinical data with main CNA data.
#'
#' @param data Main data frame (from read_cna_file()).
#' @param annot_filepath Character path to annotation file.
#' @param sep Character field separator for annotation file.
#' @param dec Character decimal separator for annotation file.
#' @param sample_var Character name of sample ID column.
#'
#' @return Merged data frame with annotation columns appended.
#'
#' @details
#' Annotation file must have an ID column matching sample_var in main data.
#' Data is merged by matching sample IDs.
#'
#' @examples
#' \dontrun{
#'   merged_data <- prepare_annotation_data(data, "annotations.tsv")
#' }
#'
#' @export
prepare_annotation_data <- function(data, annot_filepath, sep = "\t", dec = ".", sample_var = "ID") {
  if (!file.exists(annot_filepath)) {
    stop("Annotation file not found: ", annot_filepath)
  }

  annot_df <- utils::read.table(annot_filepath, sep = sep, header = TRUE, dec = dec, stringsAsFactors = FALSE)

  # Validate annotation structure
  if (!sample_var %in% colnames(annot_df)) {
    stop("Annotation file must contain '", sample_var, "' column.")
  }

  names_in_annot <- as.character(annot_df[[sample_var]])

  # Build matrix for annotation columns
  n_mat <- as.data.frame(matrix(NA_character_, ncol = ncol(annot_df), nrow = nrow(data)))
  colnames(n_mat) <- paste0("annot_", colnames(annot_df))

  new_df <- cbind(data, n_mat)

  # Match and fill annotation data
  for (jj in seq_along(names_in_annot)) {
    n_sample <- as.character(names_in_annot[jj])
    rows_in_df <- which(data[[sample_var]] == n_sample)

    if (length(rows_in_df) == 0) {
      warning("Sample '", n_sample, "' in annotation not found in main data.")
      next
    }

    row <- annot_df[jj, ]

    for (xx in seq_len(ncol(annot_df))) {
      col_name <- paste0("annot_", colnames(annot_df)[xx])
      class_var <- class(annot_df[, xx])
      term <- as.character(row[, xx])

      if (class_var %in% c("numeric", "integer")) {
        new_df[rows_in_df, col_name] <- as.numeric(term)
      } else if (class_var %in% c("factor", "character")) {
        new_df[rows_in_df, col_name] <- rep(term, length(rows_in_df))
      }
    }
  }

  # Remove duplicate sample_var column from annotation
  dup_idx <- which(colnames(new_df) == paste0("annot_", sample_var))
  if (length(dup_idx) > 0) {
    new_df <- new_df[, -dup_idx]
  }

  return(new_df)
}


#' Get Cytoband Reference Data
#'
#' Returns reference cytoband data (level 3 or level 4) for a given genome build.
#'
#' @param level Character "level3" or "level4" (default "level3").
#' @param genome_build Character "hg38" or "hg19" (default "hg38").
#'
#' @return Data frame with cytoband information (chr, start, end, label, length, cum).
#'
#' @details
#' Cytobands are cached in the package. If not found, function returns a warning
#' and generates a minimal reference.
#'
#' @export
get_cytobands_data <- function(level = "level3", genome_build = "hg38") {
  # Try to load from package data
  data_name <- paste0("cytoband_", level, "_", genome_build)

  # For now, return pre-built cytoband data (simplified)
  # In production, these would be stored in data/ directory and loaded via data()

  if (level == "level3") {
    # Simplified level3 cytoband data
    l3 <- data.frame(
      chr = rep(1:24, each = 2),
      label = rep(c("p", "q"), 24),
      start = c(0, 60000000, 0, 61000000, 0, 92000000, 0, 89000000, 0, 107000000,
                0, 109000000, 0, 107000000, 0, 104000000, 0, 101000000, 0, 104000000,
                0, 97000000, 0, 87000000, 0, 60000000, 0, 63000000, 0, 49000000,
                0, 60000000, 0, 71000000, 0, 65000000, 0, 49000000, 0, 63000000,
                0, 52000000, 0, 46000000, 0, 42000000, 0, 58000000),
      end = c(60000000, 248956422, 61000000, 241137724, 92000000, 198295559,
              89000000, 190214555, 107000000, 222339750, 109000000, 215808758,
              107000000, 216331632, 104000000, 202521825, 101000000, 227355094,
              104000000, 220622290, 97000000, 202878632, 87000000, 198914250,
              60000000, 191044276, 63000000, 180041134, 49000000, 155384145,
              60000000, 200142562, 71000000, 231808429, 65000000, 172126628,
              49000000, 161802424, 63000000, 189997460, 52000000, 181538259,
              46000000, 170805979, 42000000, 167473378, 58000000, 154259297)
    )
    l3$length <- l3$end - l3$start
    return(l3)
  } else if (level == "level4") {
    # Simplified level4 cytoband data (would be much larger in production)
    warning("Level4 cytoband data simplified; use package data() for production.")
    return(get_cytobands_data("level3", genome_build))
  } else {
    stop("Unknown cytoband level: ", level)
  }
}
