.onLoad <- function(libname, pkgname) {
  loaded_in_fork_child(detect_parallel_fork())
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("fraq ", utils::packageVersion("fraq"))
}

# We need to test if we are in a forked process, since TBB is not supported under fork
# https://www.intel.com/content/www/us/en/docs/onetbb/developer-guide-api-reference/2021-11/known-limitations.html
# Adapted from parallelly::isForkedChild():
# https://github.com/futureverse/parallelly/blob/develop/R/isForkedChild.R
detect_parallel_fork <- function() {
  if (!"parallel" %in% loadedNamespaces()) {
    return(FALSE)
  }

  ns <- asNamespace("parallel")
  is_child <- get0("isChild", mode = "function", envir = ns, inherits = FALSE)
  if (!is.function(is_child)) {
    return(FALSE)
  }

  isTRUE(is_child())
}

#' fraq: fastq processing with TBB flow graphs
#'
#' @description
#' fraq is a high-throughput toolkit for FASTQ processing. It uses a TBB flow
#' graph to coordinate concurrent I/O and compute, provides prebuilt kernels
#' for common tasks, and lets you extend it by writing custom kernels.
#'
#' @details
#' See individual function help pages for details on each kernel. The README
#' contains a short walkthrough with synthetic data.
#'
#' @keywords internal
"_PACKAGE"
