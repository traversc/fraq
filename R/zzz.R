.onAttach <- function(libname, pkgname) {
    packageStartupMessage("fraq ", utils::packageVersion("fraq"))
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
