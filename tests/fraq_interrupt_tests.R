# NOTE: run this interactively and interrupt manually with Esc or Ctrl+C\n"

if(interactive()) {
    suppressPackageStartupMessages(library(fraq))
    invisible(fraq_options("blocksize", 10L))
    input <- tempfile(fileext = ".fastq")
    output <- tempfile(fileext = ".fastq.gz")
    generate_random_fastq(
        input,
        n_reads = 2000,
        read_length = 100,
        name_prefix = "interrupt_"
    )

    cat("Interrupt test: single-threaded\n")
    fraq:::rcpp_fraq_wait(input, output, sleep_ms = 10L, nthreads = 1L)

    generate_random_fastq(
        input,
        n_reads = 2000,
        read_length = 100,
        name_prefix = "interrupt_"
    )
    cat("Interrupt test: multi-threaded\n")
    fraq:::rcpp_fraq_wait(input, output, sleep_ms = 10L, nthreads = 2L)
}