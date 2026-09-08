# Small helper shared by the original clinical/single-cell scripts.
configure_study <- function(module) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript analysis_script.R /path/to/private/settings.R")
  config_file <- normalizePath(args[1], mustWork = TRUE)
  config <- source(config_file, local = new.env(parent = globalenv()))$value
  values <- config[[module]]
  if (!is.list(values)) stop("Missing settings section: ", module)
  directory <- dirname(config_file)
  setting <- function(key) {
    value <- values[[key]]
    if (is.null(value)) stop("Configure ", module, "$", key, " in ", config_file)
    value
  }
  path <- function(key, optional = FALSE) {
    if (optional && is.null(values[[key]])) return(NA_character_)
    value <- setting(key)
    if (!is.character(value) || length(value) != 1L) stop("Expected one path for ", key)
    if (!grepl("^(/|~|[A-Za-z]:)", value)) value <- file.path(directory, value)
    path.expand(value)
  }
  assign("study_setting", setting, envir = .GlobalEnv)
  assign("study_path", path, envir = .GlobalEnv)
  set.seed(1)
  invisible(NULL)
}
