suppressPackageStartupMessages(library(fraq))

default_blocksize <- 65535L
extended <- identical(Sys.getenv("FRAQ_EXTENDED_TESTS"), "1")

if (extended) {
    n_reads <- default_blocksize * 10L
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

# fraq_concat to a FIFO whose reader attaches late. The writer cannot open the
# pipe until then, so everything is buffered in memory and only reaches the
# reader if flush() drains before closing. The input is deliberately small and
# independent of n_reads so that writing finishes well inside the reader's
# delay; otherwise the reader attaches mid-write, drains as it goes, and
# nothing is left pending at close.
concat_delay <- 2L
concat_src <- tempfile(fileext = ".fastq")
generate_random_fastq(concat_src, 5000L, 60L, "fifo_concat_")
concat_hash <- tools::md5sum(concat_src)

concat_fifo <- tempfile(fileext = ".fastq.fifo")
stopifnot(system2("mkfifo", concat_fifo) == 0L)
on.exit(unlink(concat_fifo), add = TRUE)

concat_sink <- tempfile(fileext = ".fastq")
late_drain <- processx::process$new(
    "sh",
    c("-c", sprintf("sleep %d; exec cat %s", concat_delay, shQuote(concat_fifo))),
    stdin = NULL,
    stdout = concat_sink,
    stderr = "|",
    supervise = TRUE,
    cleanup_tree = TRUE
)
on.exit(
    if (late_drain$is_alive()) late_drain$kill_tree(),
    add = TRUE,
    after = FALSE
)

fraq_concat(concat_src, concat_fifo, nthreads = threads)
late_drain$wait(timeout = 200000)
stopifnot(late_drain$get_exit_status() == 0L)
stopifnot(identical(unname(tools::md5sum(concat_sink)), unname(concat_hash)))
cat("FIFO concat with a late reader validated\n")
