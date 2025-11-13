#' Generate a random FASTQ file (optionally gzipped)
#'
#' @description
#' Creates a synthetic FASTQ file with random DNA sequences and Illumina-like
#' Phred+33 quality strings (high early-cycle quality with a gentle tail drop).
#' If `output_file` ends with `.gz`, the file is written
#' gzip-compressed via a connection.
#'
#' @details
#' Each read comprises four lines: header, sequence, `+`, and quality.
#' Headers are generated as `@<name_prefix><index>`. Sequences are sampled uniformly
#' from `ACGT`. Qualities follow a tapered profile that starts near Q37 and falls
#' toward the low 30s, with occasional low-quality spikes to mimic typical
#' Illumina output.
#'
#' @param output_file Character vector of length 1 (single-end) or 2 (paired-end)
#'   outputs. `.gz` suffixes create gzip-compressed files; otherwise plain-text
#'   FASTQ is written.
#' @param n_reads Integer number of reads to generate. Default `1e5`.
#' @param read_length Integer read length (number of bases per read). Default `100`.
#' @param name_prefix Character prefix for read names. Default `"read_"`.
#' @return Invisibly returns the path(s) written in `output_file`.
#'
#' @examples
#' # Example: small test files
#' tmp_fastq <- tempfile(fileext = ".fastq")
#' tmp_fastq_gz <- tempfile(fileext = ".fastq.gz")
#'
#' # Create plain FASTQ (500 reads, length 100)
#' generate_random_fastq(tmp_fastq, n_reads = 500, read_length = 100)
#'
#' # Create gzipped FASTQ (500 reads, length 100)
#' generate_random_fastq(tmp_fastq_gz, n_reads = 500, read_length = 100)
#'
#' # Paired-end example with overlapping mates
#' generate_random_fastq(c(tmp_fastq, tmp_fastq_gz),
#'                       n_reads = 100,
#'                       read_length = 150)
#' @importFrom stringfish random_strings
#' @export
generate_random_fastq <- function(
    output_file,
    n_reads = 1e5,
    read_length = 100,
    name_prefix = "read_"
) {
    if (missing(output_file)) {
        stop("`output_file` must be supplied.")
    }
    n_reads <- as.integer(n_reads)
    read_length <- as.integer(read_length)
    if (n_reads < 0L) {
        stop("`n_reads` must be greater than or equal to zero.")
    }
    if (read_length <= 0L) {
        stop("`read_length` must be a positive integer.")
    }
    out_paths <- as.character(output_file)
    if (!(length(out_paths) %in% c(1L, 2L))) {
        stop("`output_file` must have length 1 (single-end) or 2 (paired-end).")
    }

    seq_chars <- "ACGT"
    if (length(out_paths) == 1L) {
        write_single_fastq(
            out_paths,
            n_reads,
            read_length,
            name_prefix,
            seq_chars
        )
    } else {
        write_paired_fastq(
            out_paths,
            n_reads,
            read_length,
            name_prefix,
            seq_chars
        )
    }
}

write_single_fastq <- function(
    path,
    n_reads,
    read_length,
    name_prefix,
    seq_chars
) {
    output <- character(n_reads * 4L)
    name_idx <- seq(1L, length(output), 4L)
    seq_idx <- seq(2L, length(output), 4L)
    divider_idx <- seq(3L, length(output), 4L)
    qual_idx <- seq(4L, length(output), 4L)

    output[name_idx] <- paste0("@", name_prefix, seq_len(n_reads))
    output[seq_idx] <- stringfish::random_strings(
        N = n_reads,
        string_size = read_length,
        charset = seq_chars,
        vector_mode = "normal"
    )
    output[divider_idx] <- "+"
    output[qual_idx] <- generate_quality_strings(n_reads, read_length)

    write_lines(path, output)
    invisible(path)
}

write_paired_fastq <- function(
    paths,
    n_reads,
    read_length,
    name_prefix,
    seq_chars
) {
    overlaps <- stats::rnorm(n_reads, mean = 10, sd = 20)
    overlaps <- as.integer(round(overlaps))
    overlaps[is.na(overlaps)] <- 10L
    overlaps <- pmax(1L, pmin(read_length, overlaps))

    r1_seq <- stringfish::random_strings(
        N = n_reads,
        string_size = read_length,
        charset = seq_chars,
        vector_mode = "normal"
    )
    r2_seq_raw <- stringfish::random_strings(
        N = n_reads,
        string_size = read_length,
        charset = seq_chars,
        vector_mode = "normal"
    )
    r1_tail <- substring(r1_seq, read_length - overlaps + 1L, read_length)
    r2_suffix <- substring(r2_seq_raw, overlaps + 1L, read_length)
    r2_template <- paste0(r1_tail, r2_suffix)
    r2_seq <- reverse_complement(r2_template)

    read_ids <- seq_len(n_reads)
    qual_r1 <- generate_quality_strings(n_reads, read_length)
    qual_r2 <- generate_quality_strings(n_reads, read_length)
    r1_lines <- build_fastq_lines(
        paste0("@", name_prefix, read_ids, "/1"),
        r1_seq,
        qual_r1
    )
    r2_lines <- build_fastq_lines(
        paste0("@", name_prefix, read_ids, "/2"),
        r2_seq,
        qual_r2
    )
    write_lines(paths[[1L]], r1_lines)
    write_lines(paths[[2L]], r2_lines)
    invisible(paths)
}

generate_quality_strings <- function(
    n_reads,
    read_length,
    high_q = 37,
    tail_q = 31,
    sd = 1.5,
    low_spike_prob = 0.02,
    low_mean = 15,
    low_sd = 3
) {
    positions <- seq_len(read_length)
    frac <- if (read_length == 1L) {
        rep(0, read_length)
    } else {
        (positions - 1) / (read_length - 1)
    }
    mean_profile <- high_q - (high_q - tail_q) * (frac^1.2)

    mean_matrix <- matrix(mean_profile,
        nrow = n_reads, ncol = read_length, byrow = TRUE
    )
    noise <- matrix(
        stats::rnorm(n_reads * read_length, mean = 0, sd = sd),
        nrow = n_reads, ncol = read_length
    )
    scores <- mean_matrix + noise
    if (low_spike_prob > 0) {
        spike_mask <- matrix(
            stats::runif(n_reads * read_length) < low_spike_prob,
            nrow = n_reads, ncol = read_length
        )
        if (any(spike_mask)) {
            scores[spike_mask] <- stats::rnorm(
                sum(spike_mask),
                mean = low_mean, sd = low_sd
            )
        }
    }
    scores <- round(pmin(pmax(scores, 2), 41))
    apply(scores, 1, function(row) {
        if (!length(row)) {
            ""
        } else {
            intToUtf8(as.integer(row) + 33L)
        }
    })
}

build_fastq_lines <- function(names, seqs, quals) {
    stopifnot(length(names) == length(seqs), length(seqs) == length(quals))
    n <- length(names)
    if (n == 0L) {
        return(character())
    }
    lines <- character(n * 4L)
    idx <- seq_len(n)
    lines[(idx - 1L) * 4L + 1L] <- names
    lines[(idx - 1L) * 4L + 2L] <- seqs
    lines[(idx - 1L) * 4L + 3L] <- "+"
    lines[(idx - 1L) * 4L + 4L] <- quals
    lines
}

reverse_complement <- function(seqs) {
    if (!length(seqs)) {
        return(seqs)
    }
    complemented <- chartr("ACGTacgt", "TGCAtgca", seqs)
    vapply(
        complemented,
        function(s) {
            if (!nzchar(s)) {
                ""
            } else {
                intToUtf8(rev(utf8ToInt(s)))
            }
        },
        character(1),
        USE.NAMES = FALSE
    )
}

write_lines <- function(path, lines) {
    if (grepl("\\.zst$", path, ignore.case = TRUE)) {
        stop(
            "zst not supported in this function. Please use .gz or uncompressed."
        )
    }
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        con <- gzfile(path, "wb")
        on.exit(close(con), add = TRUE)
        writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
    } else {
        con <- file(path, open = "wb")
        on.exit(close(con), add = TRUE)
        writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
    }
}
