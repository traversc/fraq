# We need to test if we are in a forked process, since TBB is not supported under fork
# https://www.intel.com/content/www/us/en/docs/onetbb/developer-guide-api-reference/2021-11/known-limitations.html
#
# RcppParallel::isProcessForkedChild() compares the current pid against the one
# recorded when RcppParallel itself was loaded, which is the pid that owns the
# TBB runtime we use. That is the signal that actually matters, and unlike the
# old parallel::isChild() probe it does not care which package did the forking.
#
# This only covers being loaded inside an existing fork child; the C++ side
# additionally compares against fraq's own load pid to catch the more common
# case of fraq being loaded first and the fork happening afterwards.
.onLoad <- function(libname, pkgname) {
  loaded_in_fork_child(isTRUE(RcppParallel::isProcessForkedChild()))
}

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
