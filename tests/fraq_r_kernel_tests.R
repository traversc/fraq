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

cat("Testing fraq_run_r single-end identity kernel...\n")
single_records <- make_records("r_", 8L)
single_in <- write_fastq(tempfile(fileext = ".fastq"), single_records)
single_out <- tempfile(fileext = ".fastq")
seen_indices <- numeric()

fraq_run_r(
    single_in,
    function(reads, index) {
        stopifnot(identical(names(reads), "read1"))
        stopifnot(identical(names(reads[[1]]), c("name", "seq", "qual")))
        seen_indices <<- c(seen_indices, index)
        stats::setNames(list(reads[[1]]), single_out)
    },
    io_threads = 2L
)

stopifnot(identical(as.integer(seen_indices), seq_along(single_records) - 1L))
stopifnot(identical(read_fastq_records(single_out), single_records))

cat("Testing fraq_run_r NULL-returning kernel...\n")
drop_indices <- numeric()
fraq_run_r(
    single_in,
    function(reads, index) {
        drop_indices <<- c(drop_indices, index)
        NULL
    },
    io_threads = 2L
)
stopifnot(identical(as.integer(drop_indices), seq_along(single_records) - 1L))

cat("Testing fraq_run_r zero limit handling...\n")
zero_called <- FALSE
fraq_run_r(
    single_in,
    function(reads, index) {
        zero_called <<- TRUE
        NULL
    },
    limit = 0L,
    io_threads = 2L
)
stopifnot(!zero_called)

cat("Testing fraq_run_r paired-end filtering and limit...\n")
paired_r1 <- make_records("r1_", 7L)
paired_r2 <- make_records("r2_", 7L)
paired_in <- c(
    write_fastq(tempfile(fileext = ".fastq"), paired_r1),
    write_fastq(tempfile(fileext = ".fastq"), paired_r2)
)
paired_out <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))

fraq_run_r(
    paired_in,
    function(reads, index) {
        keep <- index %% 2 == 0
        stats::setNames(
            list(reads[[1]][keep, , drop = FALSE], reads[[2]][keep, , drop = FALSE]),
            paired_out
        )
    },
    limit = 5L,
    io_threads = 2L
)

expected_keep <- c(1L, 3L, 5L)
stopifnot(identical(read_fastq_records(paired_out[1]), paired_r1[expected_keep]))
stopifnot(identical(read_fastq_records(paired_out[2]), paired_r2[expected_keep]))

cat("Testing fraq_run_r demultiplexing to multiple outputs...\n")
demux_in <- write_fastq(tempfile(fileext = ".fastq"), single_records)
even_out <- tempfile(fileext = ".fastq")
odd_out <- tempfile(fileext = ".fastq")

fraq_run_r(
    demux_in,
    function(reads, index) {
        stats::setNames(
            list(
                reads[[1]][index %% 2 == 0, , drop = FALSE],
                reads[[1]][index %% 2 == 1, , drop = FALSE]
            ),
            c(even_out, odd_out)
        )
    },
    io_threads = 2L
)

stopifnot(identical(read_fastq_records(even_out), single_records[c(1L, 3L, 5L, 7L)]))
stopifnot(identical(read_fastq_records(odd_out), single_records[c(2L, 4L, 6L, 8L)]))

cat("Testing fraq_run_r .mem output...\n")
mem_key <- tempfile(fileext = ".mem")
mem_back <- tempfile(fileext = ".fastq")
fraq_run_r(
    single_in,
    function(reads, index) {
        stats::setNames(list(reads[[1]]), mem_key)
    },
    io_threads = 2L
)
fraq_convert(mem_key, mem_back, nthreads = 2L)
stopifnot(identical(read_fastq_records(mem_back), single_records))
invisible(fraq_mem_remove(mem_key))

cat("Testing fraq_run_r kernel return validation...\n")
expect_error_message(
    fraq_run_r(
        single_in,
        function(reads, index) {
            list(reads[[1]])
        },
        io_threads = 2L
    ),
    "named list"
)

cat("fraq_run_r tests completed successfully\n")
