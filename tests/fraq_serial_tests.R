suppressPackageStartupMessages(library(fraq))
library(parallel)

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

old_blocksize <- fraq_options("blocksize", 3L)
on.exit(fraq_options("blocksize", old_blocksize), add = TRUE)

cat("Testing serial C++ path with .fraq conversion...\n")
serial_in <- tempfile(fileext = ".fastq")
generate_random_fastq(
    serial_in,
    n_reads = 8L,
    read_length = 50L,
    name_prefix = "serial_"
)
serial_records <- read_fastq_records(serial_in)
serial_fraq <- tempfile(fileext = ".fraq")
serial_back <- tempfile(fileext = ".fastq")
fraq_convert(serial_in, serial_fraq, nthreads = 1L)
fraq_convert(serial_fraq, serial_back, nthreads = 1L)
stopifnot(identical(read_fastq_records(serial_back), serial_records))

cat("Testing serial concat path...\n")
concat_in <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
generate_random_fastq(
    concat_in[1],
    n_reads = 4L,
    read_length = 50L,
    name_prefix = "concat_a_"
)
generate_random_fastq(
    concat_in[2],
    n_reads = 5L,
    read_length = 50L,
    name_prefix = "concat_b_"
)
concat_a <- read_fastq_records(concat_in[1])
concat_b <- read_fastq_records(concat_in[2])
concat_out <- tempfile(fileext = ".fastq")
fraq_concat(concat_in, concat_out, nthreads = 1L)
stopifnot(identical(read_fastq_records(concat_out), c(concat_a, concat_b)))

cat("Testing serial R kernel path...\n")
paired_in <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
generate_random_fastq(
    paired_in,
    n_reads = 7L,
    read_length = 50L,
    name_prefix = "serial_r_"
)
paired_r1 <- read_fastq_records(paired_in[1])
paired_r2 <- read_fastq_records(paired_in[2])
paired_out <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
fraq_run_r(
    paired_in,
    function(reads, index) {
        keep <- index %% 2 == 0
        output <- list()
        output[[paired_out[1]]] <- reads$read1[keep, , drop = FALSE]
        output[[paired_out[2]]] <- reads$read2[keep, , drop = FALSE]
        output
    },
    nthreads = 1L
)
expected_keep <- c(1L, 3L, 5L, 7L)
stopifnot(identical(read_fastq_records(paired_out[1]), paired_r1[expected_keep]))
stopifnot(identical(read_fastq_records(paired_out[2]), paired_r2[expected_keep]))

if (.Platform$OS.type != "windows") {
    cat("Testing forked calls force serial path...\n")
    fork_results <- mclapply(seq_len(2L), function(i) {
        input <- tempfile(fileext = ".fastq")
        generate_random_fastq(
            input,
            n_reads = 6L,
            read_length = 50L,
            name_prefix = paste0("fork_", i, "_")
        )
        records <- read_fastq_records(input)
        output <- tempfile(fileext = ".fastq")
        fraq_convert(input, output, nthreads = 2L)
        convert_ok <- identical(read_fastq_records(output), records)

        r_output <- tempfile(fileext = ".fastq")
        fraq_run_r(
            input,
            function(reads, index) {
                out <- list()
                out[[r_output]] <- reads$read1
                out
            },
            nthreads = 2L
        )
        convert_ok && identical(read_fastq_records(r_output), records)
    }, mc.cores = 2L)
    stopifnot(all(unlist(fork_results, use.names = FALSE)))
} else {
    cat("Skipping forked serial guard test on this platform.\n")
}

cat("fraq serial path tests completed successfully\n")
