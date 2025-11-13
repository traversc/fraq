#' Generate an example fraq Rcpp script
#'
#' @description
#' Writes a minimal Rcpp example file showing how to write custom kernels via
#' `Rcpp`.
#'
#' @param output_file Character path where the C++ source file will be written.
#'
#' @return `NULL` invisibly.
#'
#' @examples
#' cpp <- tempfile(fileext = ".cpp")
#' fraq_rcpp_template(cpp)
#' # Rcpp::sourceCpp(cpp)  # optionally compile the example
#'
#' @export
fraq_rcpp_template <- function(output_file) {
    write(.fraq_rcpp_template_cpp, file = output_file)
}

.fraq_rcpp_template_cpp <- paste(
    c(
        "// [[Rcpp::depends(fraq)]]",
        "",
        "#include <Rcpp.h>",
        "#include <fraq.h>",
        "",
        "using fraq::Read;          // struct of strings {name, seq, qual}",
        "using fraq::output_t;      // alias for vector<pair<string, Read>>",
        "using fraq::input_t;       // alias for vector<Read>",
        "",
        "// This is a simple fraq template using Rcpp.",
        "// In R it generates a paired-end FASTQ sample, then",
        "// the kernel echoes reads to zstd-compressed outputs.",
        "// Modify process_kernel to perform your own per-record work.",
        "// You may also provide a function or functor instead of a lambda.",
        "// TBB flow handles IO and processing concurrency under the hood.",
        "",
        "// The kernel input is a vector of reads with one entry per file.",
        "// The output is a vector of {output_filename, Read} pairs.",
        "// Output size can differ from input, and file extensions",
        "// control compression for both reads and writes.",
        "",
        "// [[Rcpp::export(rng=false)]]",
        "void fraq_example(std::string input_r1, std::string input_r2,",
        "                  std::string output_r1, std::string output_r2) {",
        "    auto process_kernel = [&](input_t reads, size_t index)",
        "        -> output_t {",
        "        return {",
        "            { output_r1, std::move(reads[0]) },",
        "            { output_r2, std::move(reads[1]) }",
        "        };",
        "    };",
        "",
        "    fraq::FraqRunConfig config;",
        "    config.zstd_compress_level = 4;",
        "    fraq::run({input_r1, input_r2}, process_kernel,",
        "              5L /* nthreads */, config);",
        "}",
        "",
        "/*** R",
        "library(fraq)",
        "input_R1 <- tempfile(fileext = \".fastq.gz\")",
        "input_R2 <- tempfile(fileext = \".fastq.gz\")",
        "generate_random_fastq(input_R1, n_reads = 1000, read_length = 150,",
        "                      name_prefix = \"read1_\")",
        "generate_random_fastq(input_R2, n_reads = 1000, read_length = 150,",
        "                      name_prefix = \"read2_\")",
        "",
        "output_R1 <- tempfile(fileext = \".fastq.zst\")",
        "output_R2 <- tempfile(fileext = \".fastq.zst\")",
        "",
        "fraq_example(input_R1, input_R2, output_R1, output_R2)",
        "*/"
    ),
    collapse = "\n"
)
