suppressPackageStartupMessages(library(fraq))

write_fastq <- function(path, records) {
    stopifnot(length(records) > 0L)
    lines <- character(length(records) * 4L)
    for (i in seq_along(records)) {
        rec <- records[[i]]
        idx <- (i - 1L) * 4L
        lines[idx + 1L] <- paste0("@", rec$name)
        lines[idx + 2L] <- rec$seq
        lines[idx + 3L] <- "+"
        lines[idx + 4L] <- rec$qual
    }
    con <- file(path, open = "wb")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
    invisible(path)
}

read_fastq_records <- function(path) {
    if (!file.exists(path)) {
        return(list())
    }
    lines <- readLines(path)
    if (length(lines) == 0L) {
        return(list())
    }
    stopifnot(length(lines) %% 4L == 0L)
    n <- length(lines) / 4L
    lapply(seq_len(n), function(i) {
        idx <- (i - 1L) * 4L
        list(
            name = substring(lines[idx + 1L], 2L),
            seq = lines[idx + 2L],
            qual = lines[idx + 4L]
        )
    })
}

test_reverse_complement <- function(seq) {
    seq <- chartr("ACGTacgt", "TGCAtgca", seq)
    intToUtf8(rev(utf8ToInt(seq)))
}

compare_fastq <- function(a, b) {
    lines_a <- if (file.exists(a)) readLines(a) else character()
    lines_b <- if (file.exists(b)) readLines(b) else character()
    stopifnot(identical(lines_a, lines_b))
}

expect_error_message <- function(expr, pattern) {
    msg <- tryCatch(
        {
            expr
            NULL
        },
        error = function(e) e$message
    )
    stopifnot(!is.null(msg))
    stopifnot(grepl(pattern, msg, fixed = TRUE))
}

default_blocksize <- 65535L
invisible(fraq_options("blocksize", default_blocksize))


cat("Testing generate_random_fastq paired-end output...\n")
paired_r1 <- tempfile(fileext = ".fastq")
paired_r2 <- tempfile(fileext = ".fastq")
paired_reads <- 32L
paired_length <- 80L
generate_random_fastq(
    c(paired_r1, paired_r2),
    n_reads = paired_reads,
    read_length = paired_length,
    name_prefix = "paired_read_"
)
r1_records <- read_fastq_records(paired_r1)
r2_records <- read_fastq_records(paired_r2)
stopifnot(length(r1_records) == paired_reads)
stopifnot(length(r2_records) == paired_reads)
stopifnot(all(endsWith(vapply(r1_records, `[[`, character(1), "name"), "/1")))
stopifnot(all(endsWith(vapply(r2_records, `[[`, character(1), "name"), "/2")))
overlap_lengths <- vapply(
    seq_len(paired_reads),
    function(i) {
        r1_seq <- r1_records[[i]]$seq
        r2_seq <- r2_records[[i]]$seq
        r2_rc <- test_reverse_complement(r2_seq)
        max_k <- 0L
        for (k in paired_length:1L) {
            if (
                substr(r1_seq, paired_length - k + 1L, paired_length) ==
                    substr(r2_rc, 1L, k)
            ) {
                max_k <- k
                break
            }
        }
        max_k
    },
    integer(1)
)
stopifnot(all(overlap_lengths >= 1L))

cat("Testing fraq_options...\n")
old_val <- fraq_options("blocksize", 32L)
stopifnot(old_val == default_blocksize)
restored <- fraq_options("blocksize", default_blocksize)
stopifnot(restored == 32L)

cat("Testing fraq_downsample...\n")
downsample_records <- list(
    list(name = "r1", seq = "ACGT", qual = "IIII"),
    list(name = "r2", seq = "TGCA", qual = "HHHH"),
    list(name = "r3", seq = "GGGG", qual = "JJJJ")
)
down_in <- write_fastq(tempfile(fileext = ".fastq"), downsample_records)
invisible(fraq_options("blocksize", 4L))
down_keep <- tempfile(fileext = ".fastq")
fraq_downsample(down_in, down_keep, amount = 1, nthreads = 1L)
compare_fastq(down_in, down_keep)
set.seed(123)
down_drop <- tempfile(fileext = ".fastq")
fraq_downsample(down_in, down_drop, amount = 0, nthreads = 1L)
drop_records <- read_fastq_records(down_drop)
stopifnot(identical(length(drop_records), 0L))
invisible(fraq_options("blocksize", default_blocksize))

cat("Testing fraq_convert...\n")
conv_gz <- tempfile(fileext = ".fastq.gz")
fraq_convert(down_in, conv_gz, nthreads = 1L)
conv_back <- tempfile(fileext = ".fastq")
fraq_convert(conv_gz, conv_back, nthreads = 1L)
compare_fastq(down_in, conv_back)

cat("Testing fraq_chunk...\n")
chunk_records <- lapply(seq_len(6), function(i) {
    list(
        name = paste0("c", i),
        seq = paste(rep("A", 6), collapse = ""),
        qual = "IIIIII"
    )
})
chunk_in <- write_fastq(tempfile(fileext = ".fastq"), chunk_records)
chunk_prefix <- tempfile("fraq_chunk")
fraq_chunk(
    chunk_in,
    chunk_prefix,
    output_suffix = "fastq",
    chunk_size = 2,
    nthreads = 1L
)
chunk_size <- 2L
chunk_total <- ceiling(length(chunk_records) / chunk_size)
expected_chunks <- sprintf(
    "%s_chunk%d.fastq",
    chunk_prefix,
    0:(chunk_total - 1L)
)
stopifnot(all(file.exists(expected_chunks)))
chunk_counts <- as.integer(unname(vapply(
    expected_chunks,
    function(f) length(read_fastq_records(f)),
    integer(1)
)))
expected_counts <- rep.int(chunk_size, chunk_total)
if (chunk_total > 0L) {
    last_count <- as.integer(
        length(chunk_records) - chunk_size * (chunk_total - 1L)
    )
    expected_counts[chunk_total] <- last_count
}
stopifnot(identical(chunk_counts, expected_counts))

cat("Testing fraq_count_barcodes...\n")
barcodes <- c("ACG", "TTT", "AAC")
count_barcodes <- function(records, allow_revcomp = FALSE) {
    df <- fraq_count_barcodes(
        write_fastq(tempfile(fileext = ".fastq"), records),
        barcodes,
        max_distance = 0L,
        allow_revcomp = allow_revcomp,
        nthreads = 1L
    )
    setNames(as.integer(df$count), as.character(df$barcode))
}
count_lookup <- function(counts, key) {
    if (key %in% names(counts)) counts[[key]] else 0L
}
counts_acg <- count_barcodes(list(list(
    name = "b1",
    seq = "ACGAAAA",
    qual = "IIIIIII"
)))
stopifnot(count_lookup(counts_acg, "ACG") == 1L)
stopifnot(count_lookup(counts_acg, "NO_MATCH") == 0L)
counts_ttt <- count_barcodes(list(list(
    name = "b2",
    seq = "TTTCCCC",
    qual = "HHHHHHH"
)))
stopifnot(count_lookup(counts_ttt, "TTT") == 1L)
stopifnot(count_lookup(counts_ttt, "NO_MATCH") == 0L)
counts_revcomp <- count_barcodes(
    list(list(name = "b3", seq = "GTTGGGG", qual = "GGGGGGG")),
    allow_revcomp = TRUE
)
stopifnot(count_lookup(counts_revcomp, "AAC") == 1L)
stopifnot(count_lookup(counts_revcomp, "NO_MATCH") == 0L)
counts_revcomp_off <- count_barcodes(
    list(list(name = "b3b", seq = "GTTGGGG", qual = "GGGGGGG")),
    allow_revcomp = FALSE
)
stopifnot(count_lookup(counts_revcomp_off, "AAC") == 0L)
stopifnot(count_lookup(counts_revcomp_off, "NO_MATCH") == 1L)
counts_nomatch <- count_barcodes(list(list(
    name = "b4",
    seq = "CCCCCCC",
    qual = "IIIIIII"
)))
stopifnot(count_lookup(counts_nomatch, "NO_MATCH") == 1L)

cat("Testing barcode present in both mates does not trigger MULTI_MATCH...\n")
paired_inputs <- c(
    write_fastq(tempfile(fileext = ".fastq"), list(
        list(name = "pair/1", seq = "ACGTTTT", qual = "IIIIIII")
    )),
    write_fastq(tempfile(fileext = ".fastq"), list(
        list(name = "pair/2", seq = "GGACGTG", qual = "IIIIIII")
    ))
)
paired_df <- fraq_count_barcodes(
    paired_inputs,
    c("ACG"),
    max_distance = 0L,
    allow_revcomp = FALSE,
    nthreads = 1L
)
paired_counts <- setNames(as.integer(paired_df$count), as.character(paired_df$barcode))
stopifnot("ACG" %in% names(paired_counts))
stopifnot(paired_counts[["ACG"]] == 1L)
stopifnot(!("MULTI_MATCH" %in% paired_df$barcode))

cat("Testing fraq_demux...\n")
demux_records <- list(
    list(name = "d1", seq = "AATTTT", qual = "IIIIII"),
    list(name = "d2", seq = "CCGGGG", qual = "IIIIII"),
    list(name = "d3", seq = "GGGGGG", qual = "IIIIII")
)
demux_in <- write_fastq(tempfile(fileext = ".fastq"), demux_records)
demux_template <- file.path(
    tempdir(),
    paste0("fraq_demux_", "{barcode}", ".fastq")
)
fraq_demux(
    demux_in,
    demux_template,
    c("AA", "CC"),
    max_distance = 0L,
    nthreads = 1L
)
demux_files <- file.path(
    tempdir(),
    c("fraq_demux_AA.fastq", "fraq_demux_CC.fastq", "fraq_demux_NO_MATCH.fastq")
)
stopifnot(all(file.exists(demux_files)))
demux_counts <- vapply(
    demux_files,
    function(f) length(read_fastq_records(f)),
    integer(1)
)
stopifnot(identical(unname(demux_counts), c(1L, 1L, 1L)))

cat("Testing fraq_quality_filter...\n")
qual_records <- list(
    list(name = "q1", seq = "AAAA", qual = "IIII"),
    list(name = "q2", seq = "CCCC", qual = "!!!!"),
    list(name = "q3", seq = "GGGG", qual = "$$$$")
)
qual_in <- write_fastq(tempfile(fileext = ".fastq"), qual_records)
qual_out <- tempfile(fileext = ".fastq")
fraq_quality_filter(
    qual_in,
    qual_out,
    min_mean_quality = 30,
    max_low_q_bases = 0L,
    low_q_threshold = 20L,
    nthreads = 1L
)
qual_pass <- read_fastq_records(qual_out)
stopifnot(length(qual_pass) == 1L)
stopifnot(qual_pass[[1]]$name == "q1")

cat("Testing fraq_merge_pairs...\n")
merge_r1_records <- list(
    list(name = "pair1", seq = "ACGTAC", qual = "IIIIII"),
    list(name = "pair2", seq = "AAAAAA", qual = "IIIIII")
)
merge_r2_records <- list(
    list(name = "pair1", seq = "GTACGT", qual = "IIIIII"),
    list(name = "pair2", seq = "CCCCCC", qual = "IIIIII")
)
merge_r1 <- write_fastq(tempfile(fileext = ".fastq"), merge_r1_records)
merge_r2 <- write_fastq(tempfile(fileext = ".fastq"), merge_r2_records)
merged_out <- tempfile(fileext = ".fastq")
unmerged_out <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
merge_stats <- fraq_merge_pairs(
    c(merge_r1, merge_r2),
    merged_out,
    unmerged_out,
    min_overlap = 3L,
    max_mismatch_rate = 0.0,
    nthreads = 1L
)
merged_records <- read_fastq_records(merged_out)
stopifnot(length(merged_records) == 1L)
stopifnot(merged_records[[1]]$seq == "ACGTAC")
unmerged_records_r1 <- read_fastq_records(unmerged_out[1])
unmerged_records_r2 <- read_fastq_records(unmerged_out[2])
stopifnot(length(unmerged_records_r1) == 1L)
stopifnot(length(unmerged_records_r2) == 1L)
stopifnot(merge_stats$merged_reads == 1)
stopifnot(merge_stats$unmerged_reads == 1)

cat("Testing fraq_trim_adapters...\n")
trim_records <- list(
    list(name = "t1", seq = "AAACCC", qual = "IIIIII"),
    list(name = "t2", seq = "GGGGGG", qual = "IIIIII"),
    list(name = "t3", seq = "AAA", qual = "III")
)
trim_in <- write_fastq(tempfile(fileext = ".fastq"), trim_records)
trim_out <- tempfile(fileext = ".fastq")
trim_df <- fraq_trim_adapters(
    trim_in,
    trim_out,
    adapters = c("AAA"),
    max_distance = 0L,
    filter_untrimmed = TRUE,
    nthreads = 1L
)
trimmed_records <- read_fastq_records(trim_out)
stopifnot(length(trimmed_records) == 2L)
stopifnot(trimmed_records[[1]]$seq == "CCC")
stopifnot(trimmed_records[[2]]$seq == "")
trim_counts <- setNames(
    as.integer(trim_df$count),
    as.character(trim_df$adapter)
)
stopifnot(trim_counts[["AAA"]] == 2L)
stopifnot(trim_counts[["NO_ADAPTER"]] == 1L)

cat("Testing fraq_summary (single-end)...\n")
summary_records <- list(
    list(name = "s1", seq = "ACGT", qual = "IIII"),
    list(name = "s2", seq = "GGGG", qual = "####")
)
summary_in <- write_fastq(tempfile(fileext = ".fastq"), summary_records)
summary_se <- fraq_summary(
    summary_in,
    phred33 = TRUE,
    min_overlap = 2L,
    max_mismatch_rate = 0.0,
    nthreads = 1L
)
stopifnot(summary_se$basic_stats_R1$total_sequences[1] == 2)
stopifnot(
    summary_se$length_distribution_R1$count[
        summary_se$length_distribution_R1$length == 4
    ] ==
        2
)

cat("Testing fraq_summary (paired-end)...\n")
summary_pe <- fraq_summary(
    c(merge_r1, merge_r2),
    phred33 = TRUE,
    min_overlap = 3L,
    max_mismatch_rate = 0.0,
    nthreads = 1L
)
stopifnot(summary_pe$basic_stats_R1$total_sequences[1] == 2)
stopifnot(summary_pe$basic_stats_R2$total_sequences[1] == 2)
insert_sizes <- summary_pe$insert_size
stopifnot(is.null(insert_sizes) || any(insert_sizes$insert_size == 6))

cat("Testing fraq_summary limit...\n")
limit_records <- list(
    list(name = "l1", seq = "ACGT", qual = "IIII"),
    list(name = "l2", seq = "TGCA", qual = "HHHH"),
    list(name = "l3", seq = "GGGG", qual = "JJJJ")
)
limit_in <- write_fastq(tempfile(fileext = ".fastq"), limit_records)
limit_summary <- fraq_summary(limit_in, limit = 2L, nthreads = 1L)
stopifnot(limit_summary$basic_stats_R1$total_sequences[1] == 2)

cat("Testing fraq_slice limit and select...\n")
slice_records <- list(
    list(name = "s1", seq = "AAAA", qual = "IIII"),
    list(name = "s2", seq = "CCCC", qual = "IIII"),
    list(name = "s3", seq = "GGGG", qual = "IIII")
)
slice_in <- write_fastq(tempfile(fileext = ".fastq"), slice_records)
slice_limit_out <- tempfile(fileext = ".fastq")
fraq_slice(slice_in, slice_limit_out, limit = 2L, nthreads = 1L)
slice_limit_read <- read_fastq_records(slice_limit_out)
stopifnot(length(slice_limit_read) == 2L)
stopifnot(identical(slice_limit_read[[1]]$name, "s1"))
stopifnot(identical(slice_limit_read[[2]]$name, "s2"))

slice_select_out <- tempfile(fileext = ".fastq")
fraq_slice(slice_in, slice_select_out, select = c(0, 2), nthreads = 1L)
slice_select_read <- read_fastq_records(slice_select_out)
stopifnot(length(slice_select_read) == 2L)
stopifnot(identical(slice_select_read[[1]]$name, "s1"))
stopifnot(identical(slice_select_read[[2]]$name, "s3"))

cat("Testing fraq_slice uneven inputs emit excess read warning...\n")
uneven_r1 <- write_fastq(tempfile(fileext = ".fastq"), list(
    list(name = "u1", seq = "AAA", qual = "III"),
    list(name = "u2", seq = "CCC", qual = "III"),
    list(name = "u3", seq = "GGG", qual = "III"),
    list(name = "u4", seq = "TTT", qual = "III")
))
uneven_r2 <- write_fastq(tempfile(fileext = ".fastq"), list(
    list(name = "u1", seq = "AAA", qual = "III"),
    list(name = "u2", seq = "CCC", qual = "III")
))
uneven_out1 <- tempfile(fileext = ".fastq")
uneven_out2 <- tempfile(fileext = ".fastq")
msgs <- character()
msg_con <- textConnection("msgs", "w")
sink(msg_con, type = "message")
fraq_slice(
    c(uneven_r1, uneven_r2),
    c(uneven_out1, uneven_out2),
    limit = 3L,
    nthreads = 1L
)
sink(type = "message")
close(msg_con)
stopifnot(any(grepl("excess reads", msgs)))
uneven_reads1 <- read_fastq_records(uneven_out1)
uneven_reads2 <- read_fastq_records(uneven_out2)
stopifnot(length(uneven_reads1) == 2L)
stopifnot(length(uneven_reads2) == 2L)

cat("fraq kernel tests completed successfully\n")
