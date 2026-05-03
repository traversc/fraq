suppressPackageStartupMessages(library(fraq))

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

extended <- identical(Sys.getenv("FRAQ_EXTENDED_TESTS"), "1")
default_blocksize <- 65535L
old_blocksize <- fraq_options("blocksize", 3L)
on.exit(fraq_options("blocksize", old_blocksize), add = TRUE)

cat("Testing fraq_run_r single-end identity kernel...\n")
single_in <- tempfile(fileext = ".fastq")
generate_random_fastq(
    single_in,
    n_reads = 8L,
    read_length = 50L,
    name_prefix = "r_kernel_single_"
)
single_records <- read_fastq_records(single_in)
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
    nthreads = 2L
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
    nthreads = 2L
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
    nthreads = 2L
)
stopifnot(!zero_called)

cat("Testing fraq_run_r paired-end filtering and limit...\n")
paired_in <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
generate_random_fastq(
    paired_in,
    n_reads = 7L,
    read_length = 50L,
    name_prefix = "r_kernel_paired_"
)
paired_r1 <- read_fastq_records(paired_in[1])
paired_r2 <- read_fastq_records(paired_in[2])
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
    nthreads = 2L
)

expected_keep <- c(1L, 3L, 5L)
stopifnot(identical(read_fastq_records(paired_out[1]), paired_r1[expected_keep]))
stopifnot(identical(read_fastq_records(paired_out[2]), paired_r2[expected_keep]))

cat("Testing fraq_run_r demultiplexing to multiple outputs...\n")
demux_in <- single_in
even_out <- tempfile(fileext = ".fastq")
odd_out <- tempfile(fileext = ".fastq")
even_out_from_kernel <- even_out
odd_out_from_kernel <- odd_out
if (.Platform$OS.type == "windows") {
    even_out_from_kernel <- chartr("/", "\\", even_out)
    odd_out_from_kernel <- chartr("/", "\\", odd_out)
}

fraq_run_r(
    demux_in,
    function(reads, index) {
        stats::setNames(
            list(
                reads[[1]][index %% 2 == 0, , drop = FALSE],
                reads[[1]][index %% 2 == 1, , drop = FALSE]
            ),
            c(even_out_from_kernel, odd_out_from_kernel)
        )
    },
    nthreads = 2L
)

stopifnot(identical(read_fastq_records(even_out), single_records[c(1L, 3L, 5L, 7L)]))
stopifnot(identical(read_fastq_records(odd_out), single_records[c(2L, 4L, 6L, 8L)]))

cat("Testing fraq_run_r .mem output...\n")
mem_key <- paste0("~/", basename(tempfile("fraq_r_kernel_", fileext = ".mem")))
mem_back <- tempfile(fileext = ".fastq")
fraq_run_r(
    single_in,
    function(reads, index) {
        stats::setNames(list(reads[[1]]), mem_key)
    },
    nthreads = 2L
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
        nthreads = 2L
    ),
    "named list"
)

cat("Testing fraq_run_r R kernel errors clean up partial output...\n")
for (nt in c(1L, 2L)) {
    output_path <- tempfile(fileext = ".fastq")
    error_message <- paste0("kernel failure nthreads=", nt)
    expect_error_message(
        fraq_run_r(
            single_in,
            function(reads, index) {
                if (length(index) > 0L && index[[1L]] >= 3) {
                    stop(error_message, call. = FALSE)
                }
                stats::setNames(list(reads[[1]]), output_path)
            },
            nthreads = nt
        ),
        error_message
    )
    stopifnot(!file.exists(output_path))
}

if (extended) {
    cat("Testing fraq_run_r paired-end large identity kernel...\n")
    invisible(fraq_options("blocksize", default_blocksize))
    large_n_reads <- 1000000L
    large_input_path <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
    large_output_path <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
    generate_random_fastq(
        large_input_path,
        n_reads = large_n_reads,
        read_length = 75L,
        name_prefix = "r_kernel_large_"
    )
    large_identity_kernel <- function(reads, index) {
        stats::setNames(reads, large_output_path)
    }
    fraq_run_r(
        large_input_path,
        large_identity_kernel,
        nthreads = 4L
    )
    stopifnot(identical(
        unname(tools::md5sum(large_output_path)),
        unname(tools::md5sum(large_input_path))
    ))
}

cat("fraq_run_r tests completed successfully\n")
