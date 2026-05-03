suppressPackageStartupMessages(library(fraq))

default_blocksize <- 65535L
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

is_fifo_really_supported <- function() {
    if (!fraq_fifo_supported()) return(FALSE)
    fifo_path <- tempfile("fraq_fifo_check_", fileext = ".fifo")
    result <- tryCatch({
        status <- system(
            sprintf("mkfifo %s && echo ok || echo fail", shQuote(fifo_path)),
            intern = TRUE
        )
        unlink(fifo_path)
        identical(status, "ok")
    }, error = function(e) {
        unlink(fifo_path)
        FALSE
    })
    isTRUE(result)
}

can_run_fifo_tests <- function() {
    if (!requireNamespace("processx", quietly = TRUE)) {
        cat("Skipping FIFO test: processx not available\n")
        return(FALSE)
    }
    if (!is_fifo_really_supported()) {
        cat("Skipping FIFO test: FIFO not supported in this environment\n")
        return(FALSE)
    }
    TRUE
}

if (!can_run_fifo_tests()) {
    quit(save = "no", status = 0)
}

fifo1 <- tempfile(fileext = ".fastq.fifo")
fifo2 <- tempfile(fileext = ".fastq.fifo")
stopifnot(system2("mkfifo", c(fifo1, fifo2)) == 0L)
on.exit(unlink(c(fifo1, fifo2)), add = TRUE)

r1_src <- tempfile(fileext = ".fastq")
r2_src <- tempfile(fileext = ".fastq")
generate_random_fastq(r1_src, n_reads, read_length, "fifo_r1_")
generate_random_fastq(r2_src, n_reads, read_length, "fifo_r2_")
r1_hash <- tools::md5sum(r1_src)
r2_hash <- tools::md5sum(r2_src)

r1_sink <- tempfile(fileext = ".fastq")
r2_sink <- tempfile(fileext = ".fastq")

start_drain <- function(fifo_path, dest_path) {
    processx::process$new(
        "cat",
        c(fifo_path),
        stdin = NULL,
        stdout = dest_path,
        stderr = "|",
        supervise = TRUE,
        cleanup_tree = TRUE
    )
}

drains <- list(
    start_drain(fifo1, r1_sink),
    start_drain(fifo2, r2_sink)
)
on.exit(lapply(drains, function(p) if (p$is_alive()) p$kill_tree()), add = TRUE, after = FALSE)

fraq_convert(c(r1_src, r2_src), c(fifo1, fifo2), nthreads = threads)
cat("Done writing to FIFO, waiting on reader\n")

for (p in drains) {
    p$wait(timeout = 200000) # 200s bail out instead of hanging forever
    stopifnot(p$get_exit_status() == 0L)
}

stopifnot(identical(unname(tools::md5sum(r1_sink)), unname(r1_hash)))
stopifnot(identical(unname(tools::md5sum(r2_sink)), unname(r2_hash)))
cat("FIFO streaming validated\n")
