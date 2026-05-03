suppressPackageStartupMessages(library(fraq))
library(parallel)

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

make_records <- function(prefix, n) {
    bases <- c("ACGT", "TGCA", "GATT", "CCGG", "TTAA", "GGCC", "AACC", "TTGG")
    lapply(seq_len(n), function(i) {
        seq <- paste0(bases[((i - 1L) %% length(bases)) + 1L], bases[((i + 2L) %% length(bases)) + 1L])
        list(
            name = paste0(prefix, i),
            seq = seq,
            qual = paste(rep.int(LETTERS[((i - 1L) %% 10L) + 9L], nchar(seq)), collapse = "")
        )
    })
}

old_blocksize <- fraq_options("blocksize", 3L)
on.exit(fraq_options("blocksize", old_blocksize), add = TRUE)

cat("Testing serial C++ path with .fraq conversion...\n")
serial_records <- make_records("serial_", 8L)
serial_in <- write_fastq(tempfile(fileext = ".fastq"), serial_records)
serial_fraq <- tempfile(fileext = ".fraq")
serial_back <- tempfile(fileext = ".fastq")
fraq_convert(serial_in, serial_fraq, nthreads = 1L)
fraq_convert(serial_fraq, serial_back, nthreads = 1L)
stopifnot(identical(read_fastq_records(serial_back), serial_records))

cat("Testing serial concat path...\n")
concat_a <- make_records("concat_a_", 4L)
concat_b <- make_records("concat_b_", 5L)
concat_in <- c(
    write_fastq(tempfile(fileext = ".fastq"), concat_a),
    write_fastq(tempfile(fileext = ".fastq"), concat_b)
)
concat_out <- tempfile(fileext = ".fastq")
fraq_concat(concat_in, concat_out, nthreads = 1L)
stopifnot(identical(read_fastq_records(concat_out), c(concat_a, concat_b)))

cat("Testing serial R kernel path...\n")
paired_r1 <- make_records("serial_r1_", 7L)
paired_r2 <- make_records("serial_r2_", 7L)
paired_in <- c(
    write_fastq(tempfile(fileext = ".fastq"), paired_r1),
    write_fastq(tempfile(fileext = ".fastq"), paired_r2)
)
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
        records <- make_records(paste0("fork_", i, "_"), 6L)
        input <- write_fastq(tempfile(fileext = ".fastq"), records)
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
