suppressPackageStartupMessages(library(fraq))

read_fastq_lines <- function(path) {
    if (length(path) != 1L) {
        stop("path must be scalar", call. = FALSE)
    }
    tmp <- tempfile(fileext = ".fastq")
    on.exit(unlink(tmp), add = TRUE)
    fraq_convert(path, tmp, nthreads = 1L)
    readLines(tmp)
}

concat_identical <- function(
    inputs,
    output,
    nthreads
) {
    fraq_concat(inputs, output, nthreads = as.integer(nthreads))
    expected <- unlist(lapply(inputs, read_fastq_lines), use.names = FALSE)
    observed <- read_fastq_lines(output)
    stopifnot(identical(observed, expected))
}

set.seed(1)
extended <- identical(Sys.getenv("FRAQ_EXTENDED_TESTS"), "1")

# Test 1: plain/gz/zst serial
plain <- tempfile(fileext = ".fastq")
gz_src <- tempfile(fileext = ".fastq")
zst_src <- tempfile(fileext = ".fastq")
generate_random_fastq(
    plain,
    n_reads = 7,
    read_length = 60,
    name_prefix = "plain_"
)
generate_random_fastq(
    gz_src,
    n_reads = 11,
    read_length = 55,
    name_prefix = "gz_"
)
generate_random_fastq(
    zst_src,
    n_reads = 9,
    read_length = 45,
    name_prefix = "zst_"
)

gz_in <- tempfile(fileext = ".fastq.gz")
zst_in <- tempfile(fileext = ".fastq.zst")
fraq_convert(gz_src, gz_in, nthreads = 1L)
fraq_convert(zst_src, zst_in, nthreads = 1L)

out_fastq <- tempfile(fileext = ".fastq")
concat_identical(
    c(plain, gz_in, zst_in),
    out_fastq,
    nthreads = 1L
)

# Test 2: fraq/gz/mem parallel
fraq_src <- tempfile(fileext = ".fastq")
gz_src2 <- tempfile(fileext = ".fastq")
mem_src <- tempfile(fileext = ".fastq")
generate_random_fastq(
    fraq_src,
    n_reads = 12,
    read_length = 50,
    name_prefix = "fraq_"
)
generate_random_fastq(
    gz_src2,
    n_reads = 8,
    read_length = 52,
    name_prefix = "gz2_"
)
generate_random_fastq(
    mem_src,
    n_reads = 5,
    read_length = 48,
    name_prefix = "mem_"
)

fraq_in <- tempfile(fileext = ".fraq")
gz_in2 <- tempfile(fileext = ".fastq.gz")
mem_in <- tempfile(fileext = ".mem")
fraq_convert(fraq_src, fraq_in, nthreads = 2L)
fraq_convert(gz_src2, gz_in2, nthreads = 2L)
fraq_convert(mem_src, mem_in, nthreads = 2L)

out_fraq <- tempfile(fileext = ".fraq")
concat_identical(
    c(fraq_in, gz_in2, mem_in),
    out_fraq,
    nthreads = 3L
)

# Test 3: mixed inputs to mem output
a_src <- tempfile(fileext = ".fastq")
b_src <- tempfile(fileext = ".fastq")
generate_random_fastq(a_src, n_reads = 10, read_length = 70, name_prefix = "a_")
generate_random_fastq(b_src, n_reads = 6, read_length = 65, name_prefix = "b_")

a_fraq <- tempfile(fileext = ".fraq")
b_zst <- tempfile(fileext = ".fastq.zst")
fraq_convert(a_src, a_fraq, nthreads = 2L)
fraq_convert(b_src, b_zst, nthreads = 2L)

out_mem <- tempfile(fileext = ".mem")
concat_identical(
    c(a_fraq, b_zst),
    out_mem,
    nthreads = 2L
)

if (extended) {
    cat("Testing fraq_concat large extended input...\n")
    large_counts <- c(333333L, 333333L, 333334L)
    large_inputs <- vapply(
        seq_along(large_counts),
        function(i) {
            path <- tempfile(fileext = ".fastq")
            generate_random_fastq(
                path,
                n_reads = large_counts[i],
                read_length = 75L,
                name_prefix = paste0("concat_large_", i, "_")
            )
            path
        },
        character(1)
    )
    large_out <- tempfile(fileext = ".fastq")
    fraq_concat(large_inputs, large_out, nthreads = 4L)
    large_summary <- fraq_summary(large_out, nthreads = 1L)
    stopifnot(large_summary$basic_stats_R1$total_sequences[1] == sum(large_counts))

    large_serial_out <- tempfile(fileext = ".fastq")
    fraq_concat(large_inputs, large_serial_out, nthreads = 1L)
    stopifnot(identical(
        unname(tools::md5sum(large_serial_out)),
        unname(tools::md5sum(large_out))
    ))
}

cat("fraq_concat tests completed successfully\n")
