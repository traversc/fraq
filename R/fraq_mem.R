#' Manage in-memory FASTQ datasets
#'
#' @description
#' `fraq_mem_list()` summarizes the `.mem` datasets currently stored in the
#' session. `fraq_mem_remove()` deletes one or more `.mem` entries, freeing the
#' associated memory. `fraq_mem_load()` is a convenience wrapper around
#' [fraq_convert()] that loads on-disk FASTQ/FRAQ inputs into the in-memory
#' store after validating that the outputs end with `.mem`.
#'
#' The in-memory store lives in the current R session. For consistent results,
#' call the helper functions when no other `fraq` jobs are actively writing to
#' the same `.mem` keys.
#'
#' @return
#' * `fraq_mem_list()` returns a data frame with columns `mem_key` and
#'   `n_reads`.
#' * `fraq_mem_remove()` returns a logical vector indicating which keys were
#' removed.
#' * `fraq_mem_load()` returns the target `.mem` keys invisibly.
#'
#' @examples
#' tmp <- tempfile(fileext = ".fastq")
#' generate_random_fastq(tmp, n_reads = 100, read_length = 75)
#' mem_path <- tempfile(fileext = ".mem")
#' fraq_mem_load(tmp, mem_path)
#' fraq_mem_list()
#' fraq_mem_remove(mem_path)
#' @export
fraq_mem_list <- function() {
    rcpp_fraq_mem_list()
}

#' @rdname fraq_mem_list
#' @param mem_key Character vector of `.mem` keys to remove.
#' @export
fraq_mem_remove <- function(mem_key) {
    mem_key <- normalizePath(mem_key, winslash = "/", mustWork = FALSE)
    vapply(mem_key, rcpp_fraq_mem_remove, logical(1))
}

#' @rdname fraq_mem_list
#' @param input Character vector of FASTQ/FRAQ paths to load into memory.
#' @param nthreads Positive integer parallelism for the load.
#' @export
fraq_mem_load <- function(input, mem_key, nthreads = 1L) {
    if (!length(input)) {
        stop("`input` must contain at least one file path.")
    }
    if (length(mem_key) != length(input)) {
        stop("`mem_key` must match the length of `input`.")
    }
    if (any(!grepl("\\.mem$", mem_key, ignore.case = FALSE))) {
        stop("All `mem_key` entries must end with '.mem'.")
    }
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    mem_key <- normalizePath(mem_key, winslash = "/", mustWork = FALSE)
    rcpp_fraq_convert(input, mem_key, nthreads)
    invisible(mem_key)
}
