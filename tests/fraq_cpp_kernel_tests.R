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

is_solaris <- function() {
    !is.na(Sys.info()[["sysname"]]) &&
        grepl("SunOS", Sys.info()[["sysname"]])
}

source_cpp_with_r_tests_path <- function(code) {
    old_r_tests <- Sys.getenv("R_TESTS")
    if (nzchar(old_r_tests)) {
        Sys.setenv(
            R_TESTS = normalizePath(
                old_r_tests,
                winslash = "/",
                mustWork = TRUE
            )
        )
        on.exit(Sys.setenv(R_TESTS = old_r_tests), add = TRUE)
    }
    Rcpp::sourceCpp(code = code)
    invisible(NULL)
}

cpp_kernel_source <- '
// [[Rcpp::depends(fraq)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <fraq.h>

// [[Rcpp::export(rng=false)]]
void fraq_cpp_even_filter(std::vector<std::string> input,
                          std::vector<std::string> output,
                          int nthreads) {
    if (input.size() != output.size()) {
        Rcpp::stop("input and output must have the same length");
    }

    auto process_kernel = [&](fraq::input_t reads, std::size_t index)
        -> fraq::output_t {
        fraq::output_t out;
        if (index % 2 != 0) {
            return out;
        }
        out.reserve(reads.size());
        for (std::size_t i = 0; i < reads.size(); ++i) {
            out.push_back({output[i], std::move(reads[i])});
        }
        return out;
    };

    fraq::FraqRunConfig config;
    config.blocksize = 3;
    fraq::run(input, process_kernel, nthreads, config);
}
'

if (is_solaris()) {
    cat("Skipping custom C++ kernel sourceCpp test on Solaris.\n")
} else {
    cat("Compiling custom C++ kernel test...\n")
    source_cpp_with_r_tests_path(cpp_kernel_source)

    cat("Testing custom C++ kernel through fraq::run...\n")
    paired_in <- c(
        tempfile(fileext = ".fastq"),
        tempfile(fileext = ".fastq")
    )
    generate_random_fastq(
        paired_in,
        n_reads = 9L,
        read_length = 50L,
        name_prefix = "cpp_kernel_"
    )
    paired_r1 <- read_fastq_records(paired_in[1])
    paired_r2 <- read_fastq_records(paired_in[2])
    expected_keep <- c(1L, 3L, 5L, 7L, 9L)

    for (nt in c(1L, 2L)) {
        paired_out <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
        fraq_cpp_even_filter(paired_in, paired_out, nt)
        stopifnot(identical(
            read_fastq_records(paired_out[1]),
            paired_r1[expected_keep]
        ))
        stopifnot(identical(
            read_fastq_records(paired_out[2]),
            paired_r2[expected_keep]
        ))
    }
}

cat("fraq custom C++ kernel tests completed successfully\n")
