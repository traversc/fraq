#' Get or set FRAQ options
#'
#' @description Get or set FRAQ options.
#'
#' @param option Character string name of the option. Valid options are: 
#' "blocksize", "fraq_compress_level", "zstd_compress_level", 
#' "gzip_compress_level".
#' @param value Optional value to set the option to; if NULL, the current 
#' value is returned.
#' @return The current option value (if input `value` is `NULL`) or 
#' previous option value.
#' @examples
#' # Get current blocksize
#' fraq_options("blocksize")
#' # # Set blocksize to 16384
#' fraq_options("blocksize", 16384)
#'
fraq_options <- function(option, value = NULL) {
    if (is.null(value)) {
        rcpp_fraq_options(option, 0, FALSE)
    } else {
        rcpp_fraq_options(option, value, TRUE)
    }
}


#' Downsample FASTQ file(s)
#'
#' @description Write deterministically downsampled FASTQ file(s) to disk. 
#' `input` and `output` must be vectors of the same length (e.g., R1/R2 pairs).
#'
#' @param input Character vector of one or more input FASTQ file paths. 
#' Vectors must be the same length as `output` (e.g., R1 and R2 pairs).
#' @param output Character vector of output FASTQ file paths, same length 
#' as `input`.
#' @param amount Numeric scalar in (0, 1); proportion of reads to retain.
#' @param nthreads Integer number of threads.
#'
#' @details
#' Downsampling is deterministic: given the same inputs, `fraq_downsample()` 
#' keeps the same records every run while matching the requested proportion 
#' as closely as possible.
#' @return Invisibly returns `NULL` after writing the downsampled outputs.
#' @examples
#' r1 <- tempfile(fileext = ".fastq")
#' r2 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R1_")
#' generate_random_fastq(r2, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R2_")
#' out <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
#' fraq_downsample(c(r1, r2), out, amount = 0.1)
#' @export
fraq_downsample <- function(input, output, amount, nthreads = 1L) {
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    rcpp_fraq_downsample(input, output, amount, as.integer(nthreads))
    invisible(NULL)
}

#' Concatenate sequencing files
#'
#' @description Concatenate one or more FRAQ/FASTQ inputs (plain, 
#' `.gz`, `.zst`, `.fraq`, or `.mem`) into a single output file.
#'
#' @param input Character vector of input paths/keys to concatenate. Mixed 
#' formats are supported.
#' @param output Character scalar giving the destination path (or `.mem` key).
#' @param nthreads Integer number of threads for reading, compression, and 
#' writing.
#' @return Invisibly returns `NULL` after writing the concatenated output.
#'
#' @examples
#' tmp_dir <- tempdir()
#' inputs <- file.path(tmp_dir, sprintf("reads_%d.fastq", 1:2))
#' lapply(inputs, generate_random_fastq, n_reads = 10, read_length = 50)
#' out <- file.path(tmp_dir, "all_reads.fastq.gz")
#' fraq_concat(inputs, out, nthreads = 1L)
#' @export
fraq_concat <- function(input, output, nthreads = 1L) {
    if (length(output) != 1L) {
        stop("output must be a single file path or .mem key")
    }
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    rcpp_fraq_concat(input, output, as.integer(nthreads))
    invisible(NULL)
}

#' Convert sequencing files between supported formats
#'
#' @description Re-encode sequencing files in any supported FRAQ/FASTQ format. 
#' Input and output vectors must be the same length.
#'
#' @param input Character vector of source files/keys.
#' @param output Character vector of destination files/keys, same length 
#' as `input`.
#' @param nthreads Integer number of threads for reading/writing.
#' @details
#' FIFO pipes (paths ending with `.fifo`) are only available on Unix-like 
#' systems;
#' on Windows they are not supported and will trigger an error.
#' @return Invisibly returns `NULL` after writing the converted outputs.
#' @examples
#' src <- tempfile(fileext = ".fastq")
#' generate_random_fastq(src, n_reads = 10, read_length = 50)
#' dest <- tempfile(fileext = ".fastq.gz")
#' fraq_convert(src, dest, nthreads = 1L)
#' @export
fraq_convert <- function(input, output, nthreads = 1L) {
    if (length(input) != length(output)) {
        stop("input and output must have the same length")
    }
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    rcpp_fraq_convert(input, output, as.integer(nthreads))
    invisible(NULL)
}

#' Run an R kernel over sequencing reads
#'
#' @description Run an R function over blocks of FASTQ records. With
#' `nthreads > 1`, fraq keeps file I/O and compression work on background TBB
#' threads; with `nthreads = 1` or after fork, it uses a serial path.
#'
#' @param input Character vector of source files/keys.
#' @param kernel Function called as `kernel(reads, index)`. `reads` is a named
#' list of data frames (`read1`, `read2`, ...) with character columns `name`,
#' `seq`, and `qual`. `index` is a numeric vector of zero-based record indices
#' for the rows in each data frame.
#' @param limit Optional non-negative whole-number scalar limiting processing
#' to the first `limit` record indices. `NULL` means no limit.
#' @param nthreads Integer number of threads. When `nthreads = 1`, or when fraq
#' is running inside a forked R process, fraq uses a serial path that does not
#' construct a TBB graph. For paired-end inputs, values above 4 usually provide
#' little additional benefit.
#'
#' @details
#' The kernel must return `NULL` or a named list of data frames. Each list name
#' is an output path or `.mem` key, and each data frame must contain character
#' columns `name`, `seq`, and `qual`.
#'
#' The R kernel runs only on the calling R thread. When `nthreads > 1`,
#' background threads are used for reading, joining, demultiplexing,
#' compression, and writing.
#'
#' Blocks are delivered to the R kernel in increasing block-index order. Within
#' each call, `index` is an increasing vector of zero-based read indices.
#'
#' Do not use `parallel::mclapply()` inside the kernel. It forks the R process
#' on Unix-like systems, and forking while fraq has active background threads can
#' deadlock or crash. Use vectorized R code inside the kernel; if serial mapping
#' is needed, use `lapply()` instead.
#'
#' @return Invisibly returns `NULL` after all outputs are written.
#' @examples
#' input_paths <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
#' output_paths <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
#' generate_random_fastq(input_paths[1], n_reads = 10, read_length = 50)
#' generate_random_fastq(input_paths[2], n_reads = 10, read_length = 50)
#'
#' even_read_kernel <- function(reads, index) {
#'     keep <- index %% 2 == 0
#'     filtered_read1 <- reads$read1[keep, , drop = FALSE]
#'     filtered_read2 <- reads$read2[keep, , drop = FALSE]
#'
#'     output <- list()
#'     output[[output_paths[1]]] <- filtered_read1
#'     output[[output_paths[2]]] <- filtered_read2
#'     output
#' }
#'
#' fraq_run_r(input_paths, even_read_kernel, nthreads = 2L)
#' @export
fraq_run_r <- function(input, kernel, limit = NULL, nthreads = 1L) {
    if (!length(input)) {
        stop("`input` must contain at least one file path or .mem key.")
    }
    if (!is.function(kernel)) {
        stop("`kernel` must be a function.")
    }
    if (!is.null(limit)) {
        if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) || !is.finite(limit)) {
            stop("`limit` must be a single numeric value.")
        }
        if (limit < 0) {
            stop("`limit` must be non-negative.")
        }
        if (limit != floor(limit)) {
            stop("`limit` must be a whole number.")
        }
        if (limit == 0) {
            return(invisible(NULL))
        }
    }
    limit_val <- if (is.null(limit)) 0 else as.numeric(limit)
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    rcpp_fraq_run_r(input, kernel, limit_val, as.integer(nthreads))
    invisible(NULL)
}

#' Slice reads by index or limit
#'
#' @description Write a subset of reads from `input` to `output`, either the
#' first `limit` reads or specific zero-based indices in `select`.
#'
#' @param input Character vector of source files/keys.
#' @param output Character vector of destination files/keys, same length as
#'   `input`.
#' @param limit Numeric scalar; keep the first `limit` reads (per record index).
#'   Defaults to `NULL`.
#' @param select Numeric vector of zero-based indices to keep. Defaults to
#'   `NULL`.
#' @param nthreads Integer number of threads for reading/writing.
#'
#' @details Exactly one of `limit` or `select` must be supplied.
#'
#' @return Invisibly returns `NULL` after writing the selected reads.
#' @examples
#' src <- tempfile(fileext = ".fastq")
#' generate_random_fastq(src, n_reads = 10, read_length = 50)
#' dest <- tempfile(fileext = ".fastq")
#' fraq_slice(src, dest, limit = 5)
#' @export
fraq_slice <- function(input, output, limit = NULL, select = NULL, nthreads = 1L) {
    if (!is.null(limit) && !is.null(select)) {
        stop("Specify only one of `limit` or `select`.")
    }
    if (is.null(limit) && is.null(select)) {
        stop("Specify `limit` or `select`.")
    }
    if (!is.null(limit)) {
        if (!is.numeric(limit) || length(limit) != 1L || is.na(limit)) {
            stop("`limit` must be a single numeric value.")
        }
        if (limit < 0) {
            stop("`limit` must be non-negative.")
        }
    }
    if (!is.null(select)) {
        if (!is.numeric(select) || any(is.na(select))) {
            stop("`select` must be a numeric vector with no missing values.")
        }
        if (any(select < 0)) {
            stop("`select` must contain non-negative indices.")
        }
    }
    limit_val <- if (is.null(limit)) 0 else as.numeric(limit)
    select_vec <- if (is.null(select)) numeric() else as.numeric(select)
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    rcpp_fraq_slice(input, output, limit_val, select_vec, as.integer(nthreads))
    invisible(NULL)
}

#' Chunk sequencing files into fixed-size batches
#'
#' @description Split input datasets into sequential chunks. Each chunk is 
#' written using `output_prefix` suffixed with `_chunk{N}` and the format 
#' indicated by `output_suffix`.
#'
#' @param input Character vector of source files/keys.
#' @param output_prefix Character vector of prefixes used when naming chunked 
#' outputs. Must be the same length as `input`.
#' @param output_suffix Character scalar describing the output format; use 
#' `"fastq"`, `"fraq"`, `"mem"`, `"gz"`, `"zst"`.
#' @param chunk_size Numeric chunk size in reads; each output chunk contains 
#' up to this many records (per input stream).
#' @param nthreads Integer number of worker threads.
#'
#' @details
#' The suffix mapping follows:
#' \itemize{
#' \item `"fastq"` -> `.fastq`
#' \item `"gz"` -> `.fastq.gz`
#' \item `"zst"` -> `.fastq.zst`
#' \item `"fraq"` -> `.fraq`
#' \item `"mem"` -> `.mem`
#' }
#' @return Invisibly returns `NULL` after writing all chunked outputs.
#'
#' @examples
#' r1 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 25, read_length = 75)
#' fraq_chunk(r1,
#'            output_prefix = tempfile("chunked_R1"),
#'            output_suffix = "fastq",
#'            chunk_size = 10,
#'            nthreads = 1L)
#' @export
fraq_chunk <- function(
    input,
    output_prefix,
    output_suffix,
    chunk_size,
    nthreads = 1L
) {
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output_prefix <- normalizePath(
        output_prefix,
        winslash = "/",
        mustWork = FALSE
    )
    rcpp_fraq_chunk(
        input,
        output_prefix,
        output_suffix,
        chunk_size,
        as.integer(nthreads)
    )
    invisible(NULL)
}

#' Count barcodes in FASTQ file(s)
#'
#' @description Count occurrences of provided short barcodes in reads, 
#' allowing up to `max_distance` mismatches. Accepts one or more FASTQ files 
#' (e.g., R1/R2).
#'
#' @param input Character vector of one or more input FASTQ file paths (e.g., 
#' R1 and R2).
#' @param barcodes Character vector of barcode sequences to count.
#' @param max_distance Integer maximum number of mismatches allowed for a match.
#' @param allow_revcomp Logical; if TRUE, also match reverse complements.
#' @param nthreads Integer number of threads.
#' @return A data frame with counts per barcode.
#' @examples
#' r1 <- tempfile(fileext = ".fastq")
#' r2 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R1_")
#' generate_random_fastq(r2, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R2_")
#' short_barcodes <- c("ACGT", "TGCA", "GTAC")
#' counts <- fraq_count_barcodes(c(r1, r2), short_barcodes, max_distance = 1L, 
#'     allow_revcomp = FALSE, nthreads = 1L)
#' counts
#' @export
fraq_count_barcodes <- function(
    input,
    barcodes,
    max_distance = 1L,
    allow_revcomp = FALSE,
    nthreads = 1L
) {
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    rcpp_fraq_count_barcodes(
        input,
        barcodes,
        max_distance,
        allow_revcomp,
        nthreads
    )
}

#' Demultiplex FASTQ file(s) by barcode prefix
#'
#' @description Write barcode-specific outputs by matching a prefix on the 
#' first read of each record. Each output path is formed by substituting 
#' `{barcode}` in the supplied format string.
#'
#' @param input Character vector of one or more input FASTQ file paths 
#' (e.g., R1 and R2).
#' @param output_format Character vector of format strings, same length 
#' as `input`. Each string must contain the literal `{barcode}` placeholder.
#' @param barcodes Character vector of barcode sequences to test as prefixes.
#' @param max_distance Integer maximum Hamming distance allowed between 
#' barcode and read prefix.
#' @param nthreads Integer number of threads.
#'
#' @details
#' When no barcode matches, the literal `NO_MATCH` is substituted in place of 
#' `{barcode}`. If multiple barcodes match the same read, `MULTI_MATCH` is used.
#'
#' @return Invisibly, `NULL`. Files are written to disk according to 
#' `output_format`.
#' @examples
#' r1 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R1_")
#' out <- tempfile("R1_", fileext = "_{barcode}.fastq")
#' barcodes <- c("ACGT", "TGCA", "GTAC")
#' fraq_demux(r1, out, barcodes, max_distance = 1L, nthreads = 1L)
#' @export
fraq_demux <- function(
    input,
    output_format,
    barcodes,
    max_distance = 1L,
    nthreads = 1L
) {
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output_format <- normalizePath(
        output_format,
        winslash = "/",
        mustWork = FALSE
    )
    rcpp_fraq_demux(
        input,
        output_format,
        barcodes,
        as.integer(max_distance),
        as.integer(nthreads)
    )
    invisible(NULL)
}

#' Filter reads by whole-read quality
#'
#' @description Drop read sets when any mate fails the quality thresholds. 
#' Qualities are interpreted as PHRED+33.
#'
#' @param input Character vector of one or more input FASTQ file paths. Must 
#' be the same length as `output`.
#' @param output Character vector of output FASTQ paths.
#' @param min_mean_quality Numeric minimum mean base quality required for 
#' each mate.
#' @param max_low_q_bases Integer maximum number of bases below 
#' `low_q_threshold` allowed per mate.
#' @param low_q_threshold Integer PHRED cutoff used to count low-quality bases.
#' @param nthreads Integer number of worker threads.
#'
#' @details
#' Both thresholds are evaluated separately on every mate. If any mate fails, 
#' the entire read set is discarded.
#'
#' @return Invisibly, `NULL`. Reads are written to `output` paths for records 
#' that pass the filters.
#' @examples
#' r1 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R1_")
#' out <- tempfile(fileext = ".fastq")
#' fraq_quality_filter(r1, out, min_mean_quality = 25, max_low_q_bases = 2L, 
#'     low_q_threshold = 20L, nthreads = 1L)
#' @export
fraq_quality_filter <- function(
    input,
    output,
    min_mean_quality = 20,
    max_low_q_bases = .Machine$integer.max,
    low_q_threshold = 20L,
    nthreads = 1L
) {
    if (is.infinite(max_low_q_bases)) {
        max_low_q_bases <- as.integer(.Machine$integer.max)
    } else {
        max_low_q_bases <- as.integer(max_low_q_bases)
    }
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    rcpp_fraq_quality_filter(
        input,
        output,
        as.numeric(min_mean_quality),
        max_low_q_bases,
        as.integer(low_q_threshold),
        as.integer(nthreads)
    )
    invisible(NULL)
}

#' Merge paired-end reads into a consensus
#'
#' @description Merge R1/R2 pairs by overlapping sequences (optionally 
#' reverse-complementing R2), emitting merged reads and optional unmerged outputs.
#'
#' @param input Character vector of length 2 with input FASTQ paths (R1, R2).
#' @param output_merged Character scalar path/key receiving merged single-end 
#' reads.
#' @param output_unmerged Optional character vector of length 2 for unmerged 
#' R1/R2 outputs. Use `NULL` to drop unmerged pairs.
#' @param min_overlap Integer minimum overlap required to attempt merging.
#' @param max_mismatch_rate Numeric maximum mismatch fraction allowed within 
#' the overlap.
#' @param consensus_mode Character string controlling consensus base 
#' selection: `"max"`, `"mean"`, `"r1"`, or `"r2"`.
#' @param trim_overhang Logical; if `TRUE`, include non-overlapping tails when 
#' constructing the merged read.
#' @param revcomp_R2 Logical; if `TRUE`, reverse-complement R2 before merging.
#' @param nthreads Integer number of worker threads.
#'
#' @details Qualities are interpreted as PHRED+33.
#'
#' @return A list summarising merge statistics (`merged_reads`, 
#'     `unmerged_reads`, `mean_insert_size`, etc.).
#' @examples
#' r1 <- tempfile(fileext = ".fastq")
#' r2 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 100, read_length = 100, 
#'     name_prefix = "read_R1_")
#' generate_random_fastq(r2, n_reads = 100, read_length = 100, 
#'     name_prefix = "read_R2_")
#' out_merged <- tempfile(fileext = ".fastq")
#' fraq_merge_pairs(c(r1, r2), out_merged, output_unmerged = NULL, 
#'     min_overlap = 20L, max_mismatch_rate = 0.05)
#' @export
fraq_merge_pairs <- function(
    input,
    output_merged,
    output_unmerged = NULL,
    min_overlap = 12L,
    max_mismatch_rate = 0.10,
    consensus_mode = c("max", "mean", "r1", "r2"),
    trim_overhang = TRUE,
    revcomp_R2 = TRUE,
    nthreads = 1L
) {
    if (!is.character(input) || length(input) != 2L) {
        stop("input must be a character vector of length 2 (R1, R2)")
    }
    if (!is.character(output_merged) || length(output_merged) != 1L) {
        stop("output_merged must be a character scalar")
    }
    if (is.null(output_unmerged)) {
        output_unmerged_vec <- character()
    } else {
        if (!is.character(output_unmerged) || length(output_unmerged) != 2L) {
            stop(
                "output_unmerged must be NULL or character vector of length 2"
            )
        }
        output_unmerged_vec <- output_unmerged
    }
    consensus_mode <- match.arg(consensus_mode)
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output_merged <- normalizePath(
        output_merged,
        winslash = "/",
        mustWork = FALSE)
    output_unmerged_vec <- normalizePath(
        output_unmerged_vec,
        winslash = "/",
        mustWork = FALSE)
    rcpp_fraq_merge_pairs(
        input,
        output_merged,
        output_unmerged_vec,
        as.integer(min_overlap),
        as.numeric(max_mismatch_rate),
        consensus_mode,
        as.logical(trim_overhang),
        as.logical(revcomp_R2),
        as.integer(nthreads))
}

#' Trim adapters from FASTQ file(s)
#'
#' @description Trim occurrences of adapter sequence(s) at the start of the 
#' first fastq input. `input` and `output` must be vectors of the same 
#' length (e.g., R1/R2 pairs).
#' Adapters will be trimmed only for the first fastq, but all inputs will be 
#' filtered if `filter_untrimmed` is `TRUE`.
#'
#' @param input Character vector of one or more input FASTQ file paths. 
#' Vectors must be the same length as `output` (e.g., R1 and R2 pairs).
#' @param output Character vector of output FASTQ file paths, same length as 
#' `input`.
#' @param adapters Character vector of adapter sequences to trim. Adapters are 
#' given priority based on the order they appear.
#' @param max_distance Integer maximum number of mismatches for adapter matching.
#' @param filter_untrimmed Logical; if TRUE, drop reads with no trim.
#' @param nthreads Integer number of threads.
#' @return A data frame of counts of trimmed adapters.
#' @examples
#' r1 <- tempfile(fileext = ".fastq")
#' r2 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R1_")
#' generate_random_fastq(r2, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R2_")
#' out <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
#' adapters <- c("ACGT", "TGCA", "GTAC")
#' fraq_trim_adapters(c(r1, r2), out, adapters, max_distance = 1L, 
#'     filter_untrimmed = TRUE, nthreads = 1L)
#' @export
fraq_trim_adapters <- function(
    input,
    output,
    adapters,
    max_distance = 1L,
    filter_untrimmed = TRUE,
    nthreads = 1L
) {
    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    rcpp_fraq_trim_adapters(
        input,
        output,
        adapters,
        max_distance,
        filter_untrimmed,
        nthreads
    )
}

#' Summarize FASTQ quality metrics (single- or paired-end)
#'
#' @description
#' Compute QC summaries for single- or paired-end FASTQ files.
#' When two inputs are provided, R1 and R2 are summarized separately and an 
#' insert-size
#' histogram is reported (estimated from R1 vs reverse-complemented R2 
#' overlap).
#'
#' @param input Character vector of length 1 or 2 with input FASTQ paths.
#'   Length 1 = single-end; length 2 = paired-end (first element maps to R1, 
#' second to R2).
#' @param phred33 Logical; TRUE if qualities are PHRED+33, FALSE for PHRED+64.
#' @param min_overlap Integer minimum overlap used for insert-size estimation 
#' (paired-end only).
#' @param max_mismatch_rate Numeric between 0 and 1 (inclusive); maximum 
#' allowed mismatch rate within the overlapped region (paired-end only).
#' @param limit Numeric cap on the number of read sets to process. Use `0` to
#'   process all available reads.
#' @param nthreads Integer number of threads.
#'
#' @details
#' Outputs per-mate tables:
#' \itemize{
#' \item \code{basic_stats_R{1,2}}: total sequences, total bases, min/mean/max 
#' length, GC percent.
#' \item \code{per_base_quality_R{1,2}}: mean PHRED by 1-based position (with 
#' counts).
#' \item \code{per_base_content_R{1,2}}: long format base usage by position 
#' (A/C/G/T/N/other).
#' \item \code{length_distribution_R{1,2}}: histogram of sequence lengths.
#' \item \code{avg_read_quality_R{1,2}}: histogram of rounded per-read average
#' quality (columns \code{avg_quality}, \code{count}).
#' }
#' For paired-end inputs, \code{insert_size} is included when overlaps are 
#' found.
#'
#' @return A named list of data frames. For single-end: R1-only tables. For 
#' paired-end: R1/R2 tables plus optional \code{insert_size}. Each mate includes
#' \code{basic_stats}, per-base quality/content, length distributions, and
#' average read quality histograms.
#'
#' @examples
#' # Single-end example
#' r1 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R1_")
#' res_se <- fraq_summary(r1, nthreads = 1L)
#' res_se$basic_stats_R1
#'
#' # Paired-end example
#' r1 <- tempfile(fileext = ".fastq"); r2 <- tempfile(fileext = ".fastq")
#' generate_random_fastq(r1, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R1_")
#' generate_random_fastq(r2, n_reads = 1000, read_length = 100, 
#'     name_prefix = "read_R2_")
#' res_pe <- fraq_summary(c(r1, r2), nthreads = 1L)
#' res_pe$basic_stats_R1
#' # dplyr example using pipes
#' # library(dplyr)
#' # res_pe$insert_size %>% arrange(desc(count)) %>% head()
#'
#' @export
fraq_summary <- function(
    input,
    phred33 = TRUE,
    min_overlap = 12L,
    max_mismatch_rate = 0.10,
    limit = 0L,
    nthreads = 1L
) {
    if (!is.character(input) || length(input) < 1L || length(input) > 2L) {
        stop("input must be a character vector of length 1 or 2")
    }
    if (!is.logical(phred33) || length(phred33) != 1L) {
        stop("phred33 must be a single logical")
    }
    if (
        !is.numeric(max_mismatch_rate) ||
            max_mismatch_rate < 0 ||
            max_mismatch_rate > 1
    ) {
        stop("max_mismatch_rate must be a number in [0,1]")
    }
    if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) || limit < 0) {
        stop("limit must be a single non-negative number")
    }

    input <- normalizePath(input, winslash = "/", mustWork = FALSE)
    rcpp_fraq_summary(
        input,
        phred33,
        as.integer(min_overlap),
        as.numeric(max_mismatch_rate),
        as.numeric(limit),
        as.integer(nthreads)
    )
}


#' Align a query to a target
#'
#' @description Calculate distances between query sequences and a target under 
#' a chosen boundary model using Levenshtein or Hamming distance.
#'
#' @param query Character vector or Biostrings `XString`/`XStringSet` of 
#' query sequences.
#' @param target Character vector or Biostrings `XString`/`XStringSet` of 
#' target sequences.
#' @param max_distance Integer maximum allowed distance; defaults to 
#' .Machine$integer.max.
#' @param ambiguity_base Single character to treat as ambiguity when 
#' matching, or empty string "" to disable; must be length 0 or 1.
#' @param boundary One of "contains", "global", or "starts".
#' @param distance_metric One of "lv" (Levenshtein) or "hm" (Hamming). "hm" 
#' requires query and target to be the same length.
#' @return A data frame with the Biostrings inputs as the first two columns 
#' followed by the alignment metadata.
#' @examples
#' fraq_align("ACGTNT", "ACGTAT", max_distance = 2L, ambiguity_base = "N",
#'            boundary = "contains", distance_metric = "lv")
#' fraq_align(Biostrings::DNAString("ACGT"), Biostrings::DNAString("ACGA"),
#'            max_distance = 1L, boundary = "global", distance_metric = "hm")
#' @importFrom Biostrings BStringSet
#' @export
fraq_align <- function(
    query,
    target,
    max_distance = 2147483647L,
    ambiguity_base = "",
    boundary = "contains",
    distance_metric = "lv"
) {
    normalize_fraq_align_input <- function(x, arg_name) {
        if (is.character(x)) {
            return(x)
        }
        if (inherits(x, "XStringSet") || inherits(x, "XString")) {
            return(as.character(x))
        }
        stop(
            sprintf(
            "`%s` must be a character vector or Biostrings XString/XStringSet.",
            arg_name
            ),
            call. = FALSE
        )
    }
    query_chr <- normalize_fraq_align_input(query, "query")
    target_chr <- normalize_fraq_align_input(target, "target")
    align_df <- rcpp_fraq_align(
        query_chr,
        target_chr,
        max_distance,
        ambiguity_base,
        boundary,
        distance_metric
    )
    query_bs <- Biostrings::BStringSet(query_chr)
    target_bs <- Biostrings::BStringSet(target_chr)
    align_df$query <- query_bs
    align_df$target <- target_bs
    cols <- c("query", "target", setdiff(names(align_df), 
        c("query", "target")))
    align_df[, cols, drop = FALSE]
}

#' Detect FRAQ FIFO support
#'
#' @description Report whether the current build of **fraq** was compiled with
#'     named pipe (FIFO) support. FIFO outputs (paths ending in `.fifo`) are only
#'     available on Unix-like platforms where the build detected `S_IFIFO`.
#'
#' @details The result is determined at compile time; reinstalling the package
#'     on an operating system that exposes FIFOs is required to enable support.
#'
#' @return Logical scalar indicating whether FIFO inputs/outputs are supported.
#'
#' @examples
#' if (fraq_fifo_supported()) {
#'     message("FIFO streams are available on this platform.")
#' } else {
#'     message("Use regular files instead of .fifo paths on this build.")
#' }
#'
#' @name fraq_fifo_supported
#' @aliases fraq_fifo_supported
#' @usage fraq_fifo_supported()
NULL
