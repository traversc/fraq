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

uint_le_raw <- function(value, bytes) {
    out <- raw(bytes)
    value <- as.double(value)
    for (i in seq_len(bytes)) {
        out[i] <- as.raw(value %% 256)
        value <- floor(value / 256)
    }
    out
}

fraq_file_header_raw <- function() {
    c(charToRaw("FRAQ"), raw(10L), as.raw(0L), as.raw(1L))
}

write_raw_fraq <- function(bytes) {
    path <- tempfile(fileext = ".fraq")
    con <- file(path, open = "wb")
    on.exit(close(con), add = TRUE)
    writeBin(bytes, con)
    path
}

write_fraq_with_minimal_count_header <- function(count_raw, count_class) {
    write_raw_fraq(c(
        fraq_file_header_raw(),
        uint_le_raw(count_class, 4L),
        count_raw,
        raw(8L)
    ))
}

default_blocksize <- 65535L
invisible(fraq_options("blocksize", default_blocksize))

extended <- identical(Sys.getenv("FRAQ_EXTENDED_TESTS"), "1")

if (extended) {
    n_reads <- 1000000L
    read_length <- 100L
    threads <- 4L
    invisible(fraq_options("blocksize", default_blocksize))
} else {
    n_reads <- 24L
    read_length <- 60L
    threads <- 1L
    small_block <- 7L
    invisible(fraq_options("blocksize", small_block))
}

expect_fraq_convert_error <- function(path, pattern) {
    for (nt in unique(c(1L, threads))) {
        expect_error_message(
            fraq_convert(path, tempfile(fileext = ".fastq"), nthreads = nt),
            pattern
        )
    }
}

cat("Generating source FASTQ for format tests...\n")
src_fastq <- tempfile(fileext = ".fastq")
generate_random_fastq(
    src_fastq,
    n_reads = n_reads,
    read_length = read_length,
    name_prefix = "fmt_"
)
hash_original <- tools::md5sum(src_fastq)

formats <- c(".fastq", ".fastq.gz", ".fastq.zst", ".fraq", ".mem")

for (fmt in formats) {
    cat("Testing round-trip conversion through format:", fmt, "\n")
    target <- tempfile(fileext = fmt)
    fraq_convert(src_fastq, target, nthreads = threads)
    back <- tempfile(fileext = ".fastq")
    fraq_convert(target, back, nthreads = threads)
    hash_back <- tools::md5sum(back)
    stopifnot(identical(as.character(hash_back), as.character(hash_original)))
}


cat("Testing format signature validation...\n")
bogus_gz <- tempfile(fileext = ".fastq.gz")
writeBin(charToRaw("not a gzip"), bogus_gz)
expect_error_message(
    fraq_convert(bogus_gz, tempfile(fileext = ".fastq"), nthreads = threads),
    "valid gzip"
)

bogus_zst <- tempfile(fileext = ".fastq.zst")
writeBin(charToRaw("not a zstd file"), bogus_zst)
expect_error_message(
    fraq_convert(bogus_zst, tempfile(fileext = ".fastq"), nthreads = threads),
    "valid zstd"
)

bogus_fraq <- tempfile(fileext = ".fraq")
writeBin(raw(4), bogus_fraq)
expect_error_message(
    fraq_convert(bogus_fraq, tempfile(fileext = ".fastq"), nthreads = threads),
    "FRAQ"
)

reserved_header <- fraq_file_header_raw()
reserved_header[5L] <- as.raw(1L)
expect_fraq_convert_error(
    write_raw_fraq(reserved_header),
    "unsupported FRAQ header flags"
)

big_endian_header <- fraq_file_header_raw()
big_endian_header[15L] <- as.raw(1L)
expect_fraq_convert_error(
    write_raw_fraq(big_endian_header),
    "big-endian FRAQ files are not supported"
)

expect_fraq_convert_error(
    write_raw_fraq(c(fraq_file_header_raw(), as.raw(c(0L, 0L, 0L, 128L)))),
    "unsupported FRAQ block metadata flags"
)

invalid_counts <- list(
    write_fraq_with_minimal_count_header(as.raw(0L), 0L),
    write_fraq_with_minimal_count_header(uint_le_raw(65536L, 4L), 2L),
    write_fraq_with_minimal_count_header(as.raw(rep(255L, 8L)), 3L)
)
for (path in invalid_counts) {
    expect_fraq_convert_error(
        path,
        "invalid FRAQ block read count"
    )
}

seq_size_src <- write_fastq(
    tempfile(fileext = ".fastq"),
    list(
        list(name = "seq_size", seq = "ACGT", qual = "IIII")
    )
)
seq_size_fraq <- tempfile(fileext = ".fraq")
fraq_convert(seq_size_src, seq_size_fraq, nthreads = 1L)
seq_size_bytes <- readBin(seq_size_fraq, what = "raw", n = file.info(seq_size_fraq)$size)
seq_size_bytes[23L] <- as.raw(1L)
expect_fraq_convert_error(
    write_raw_fraq(seq_size_bytes),
    "sequence payload size mismatch"
)

qual_mismatch <- write_fastq(
    tempfile(fileext = ".fastq"),
    list(
        list(name = "qual_mismatch", seq = "ACGT", qual = "IIIIII")
    )
)
for (nt in unique(c(1L, threads))) {
    output_path <- tempfile(fileext = ".fraq")
    expect_error_message(
        fraq_convert(qual_mismatch, output_path, nthreads = nt),
        "sequence and quality lengths differ"
    )
    stopifnot(!file.exists(output_path))
}
for (nt in unique(c(1L, threads))) {
    output_path <- tempfile(fileext = ".fraq")
    expect_error_message(
        fraq_concat(qual_mismatch, output_path, nthreads = nt),
        "sequence and quality lengths differ"
    )
    stopifnot(!file.exists(output_path))
}

cat("Testing fraq_mem helpers...\n")
mem_src <- tempfile(fileext = ".fastq")
generate_random_fastq(
    mem_src,
    n_reads = 12,
    read_length = 40,
    name_prefix = "mem_src_"
)
mem_key <- tempfile(fileext = ".mem")
mem_key <- normalizePath(mem_key, winslash = "/", mustWork = FALSE)
fraq_mem_load(mem_src, mem_key, nthreads = 1L)
mem_df <- fraq_mem_list()
stopifnot(mem_key %in% mem_df$mem_key)
row_idx <- which(mem_df$mem_key == mem_key)
stopifnot(length(row_idx) == 1L)
stopifnot(as.integer(mem_df$n_reads[row_idx]) == 12L)
mem_removed <- fraq_mem_remove(c(mem_key, "nonexistent.mem"))
stopifnot(identical(unname(mem_removed), c(TRUE, FALSE)))
mem_df2 <- fraq_mem_list()
stopifnot(!(mem_key %in% mem_df2$mem_key))

cat("Testing fraq_mem replacement and multi-load...\n")
mem_src_new <- tempfile(fileext = ".fastq")
generate_random_fastq(
    mem_src_new,
    n_reads = 7,
    read_length = 35,
    name_prefix = "mem_replace_"
)
fraq_mem_load(mem_src_new, mem_key, nthreads = 1L)
mem_df3 <- fraq_mem_list()
new_count <- as.integer(mem_df3$n_reads[mem_df3$mem_key == mem_key])
stopifnot(new_count == 7L)

multi_inputs <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
generate_random_fastq(
    multi_inputs[1],
    n_reads = 3,
    read_length = 20,
    name_prefix = "multi1_"
)
generate_random_fastq(
    multi_inputs[2],
    n_reads = 3,
    read_length = 25,
    name_prefix = "multi2_"
)
multi_keys <- c(tempfile(fileext = ".mem"), tempfile(fileext = ".mem"))
multi_keys <- normalizePath(multi_keys, winslash = "/", mustWork = FALSE)
fraq_mem_load(multi_inputs, multi_keys, nthreads = 1L)
mem_df4 <- fraq_mem_list()
stopifnot(all(mapply(
    function(k, n) {
        entries <- mem_df4[mem_df4$mem_key == k, , drop = FALSE]
        length(entries$n_reads) == 1L && as.integer(entries$n_reads) == n
    },
    multi_keys,
    c(3L, 3L)
)))

cat("Testing argument validation...\n")
invalid_in <- write_fastq(
    tempfile(fileext = ".fastq"),
    list(
        list(name = "r1", seq = "ACGT", qual = "IIII")
    )
)
invalid_out <- tempfile(fileext = ".fastq")
expect_error_message(
    fraq_downsample(invalid_in, invalid_out, amount = 1.5, nthreads = 1L),
    "amount must between 0 and 1"
)

concat_a <- tempfile(fileext = ".fastq")
concat_b <- tempfile(fileext = ".fastq")
generate_random_fastq(
    concat_a,
    n_reads = 3,
    read_length = 20,
    name_prefix = "ca_"
)
generate_random_fastq(
    concat_b,
    n_reads = 2,
    read_length = 20,
    name_prefix = "cb_"
)
expect_error_message(
    fraq_concat(
        c(concat_a, concat_b),
        c("out1.fastq", "out2.fastq"),
        nthreads = 1L
    ),
    "output must be a single file path or .mem key"
)

mem_len_mismatch <- tempfile(fileext = ".fastq")
generate_random_fastq(
    mem_len_mismatch,
    n_reads = 2,
    read_length = 20,
    name_prefix = "mlm_"
)
expect_error_message(
    fraq_mem_load(
        c(mem_len_mismatch, mem_len_mismatch),
        tempfile(fileext = ".mem"),
        nthreads = 1L
    ),
    "`mem_key` must match the length of `input`."
)

cat("Testing reader imbalance warnings...\n")
imbal_r1 <- tempfile(fileext = ".fastq")
imbal_r2 <- tempfile(fileext = ".fastq")
generate_random_fastq(
    imbal_r1,
    n_reads = 6,
    read_length = 40,
    name_prefix = "imb1_"
)
generate_random_fastq(
    imbal_r2,
    n_reads = 3,
    read_length = 40,
    name_prefix = "imb2_"
)
imbal_out1 <- tempfile(fileext = ".fastq")
imbal_out2 <- tempfile(fileext = ".fastq")
imbal_msg <- capture.output(
    fraq_convert(
        c(imbal_r1, imbal_r2),
        c(imbal_out1, imbal_out2),
        nthreads = 1L
    ),
    type = "message"
)
stopifnot(any(grepl("inputs: 0", imbal_msg)))

balanced_r1 <- tempfile(fileext = ".fastq")
balanced_r2 <- tempfile(fileext = ".fastq")
generate_random_fastq(
    balanced_r1,
    n_reads = 4,
    read_length = 30,
    name_prefix = "bal1_"
)
generate_random_fastq(
    balanced_r2,
    n_reads = 4,
    read_length = 30,
    name_prefix = "bal2_"
)
balanced_out1 <- tempfile(fileext = ".fastq")
balanced_out2 <- tempfile(fileext = ".fastq")
bal_msg <- capture.output(
    fraq_convert(
        c(balanced_r1, balanced_r2),
        c(balanced_out1, balanced_out2),
        nthreads = 1L
    ),
    type = "message"
)
stopifnot(length(bal_msg) == 0L)

cat("fraq format round-trip tests completed successfully\n")
