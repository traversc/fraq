#' Bridge FRAQ formats with ShortReadQ
#'
#' @description
#' `fraq_export_shortreadq()` converts FRAQ/FASTQ inputs into in-memory
#' `ShortReadQ` objects via `fraq_convert`. `fraq_import_shortreadq`
#' converts `ShortReadQ` objects to any supported FRAQ format.
#'
#' @param input Character vector of input paths/keys accepted by
#'   [fraq_convert()].
#' @param nthreads Positive integer passed to `fraq_convert`.
#' @param tmpdir Directory used for staging temporary files.
#' @param shortreadq A `ShortReadQ` or list of `ShortReadQ` objects to be 
#' written via FRAQ encoders.
#' @param output Character vector of destination paths/keys, same length as
#'   `shortreadq`.
#'
#' @return
#' * `fraq_export_shortreadq()` returns a single `ShortReadQ` when `input` has
#'   length 1, otherwise a list of `ShortReadQ` objects.
#' * `fraq_import_shortreadq()` invisibly returns the normalized `output`
#'   vector after conversion.
#'
#' @examples
#' if (requireNamespace("ShortRead", quietly = TRUE)) {
#' fq <- tempfile(fileext = ".fastq")
#' generate_random_fastq(fq, n_reads = 10, read_length = 50)
#' fraq_path <- tempfile(fileext = ".fraq")
#' fraq_convert(fq, fraq_path)
#'
#' reads <- fraq_export_shortreadq(fraq_path)
#' roundtrip_fastq <- tempfile(fileext = ".fastq")
#' fraq_import_shortreadq(reads, roundtrip_fastq)
#' stopifnot(file.exists(roundtrip_fastq))
#' }
#'
#' @seealso `ShortRead::readFastq()`, `ShortRead::writeFastq()`,
#'   [fraq_convert()]
#' @name fraq_shortread
NULL

check_shortread_installed <- function(caller) {
    if (!requireNamespace("ShortRead", quietly = TRUE)) {
        stop(
            sprintf(
                "`%s()` requires the suggested package ShortRead.",
                caller
            ),
            "\nInstall it with BiocManager::install(\"ShortRead\").",
            call. = FALSE
        )
    }
}

#' @rdname fraq_shortread
#' @export
fraq_export_shortreadq <- function(input, nthreads = 1L, tmpdir = tempdir()) {
    check_shortread_installed("fraq_export_shortreadq")
    if (!length(input)) {
        stop("`input` must contain at least one path or key.")
    }
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    tmp_paths <- vapply(
        seq_along(input),
        function(i) {
            tempfile(
                sprintf("fraq_export_%d_", i),
                tmpdir = tmpdir,
                fileext = ".fastq"
            )
        },
        character(1),
        USE.NAMES = FALSE
    )
    fraq_convert(input, tmp_paths, nthreads = nthreads)
    results <- lapply(tmp_paths, ShortRead::readFastq)
    unlink(tmp_paths, recursive = FALSE, force = FALSE)
    if (length(results) == 1L) {
        results[[1L]]
    } else {
        results
    }
}

#' @rdname fraq_shortread
#' @export
fraq_import_shortreadq <- function(
    shortreadq,
    output,
    nthreads = 1L,
    tmpdir = tempdir()
) {
    check_shortread_installed("fraq_import_shortreadq")
    if (inherits(shortreadq, "ShortReadQ")) {
        sr_list <- list(shortreadq)
    } else if (is.list(shortreadq) && length(shortreadq)) {
        ok <- vapply(shortreadq, inherits, logical(1), "ShortReadQ")
        if (!all(ok)) {
            stop(
                "All elements of `shortreadq` must be ShortReadQ objects.",
                call. = FALSE
            )
        }
        sr_list <- shortreadq
    } else {
        stop(
            "`shortreadq` must be a ShortReadQ or list of ShortReadQ objects.",
            call. = FALSE
        )
    }
    if (length(sr_list) != length(output)) {
        stop("`shortreadq` and `output` must have the same length.")
    }
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    tmp_paths <- vapply(
        seq_along(sr_list),
        function(i) {
            tempfile(
                sprintf("fraq_import_%d_", i),
                tmpdir = tmpdir,
                fileext = ".fastq"
            )
        },
        character(1),
        USE.NAMES = FALSE
    )
    for (i in seq_along(sr_list)) {
        ShortRead::writeFastq(sr_list[[i]], tmp_paths[[i]], compress = FALSE)
    }
    fraq_convert(tmp_paths, output, nthreads = nthreads)
    unlink(tmp_paths, recursive = FALSE, force = FALSE)
    invisible(output)
}
