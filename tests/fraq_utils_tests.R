suppressPackageStartupMessages(library(fraq))
suppressPackageStartupMessages(library(Biostrings))

compare_fastq_records <- function(src_path, converted_path, nthreads = 1L) {
    src_tmp <- tempfile(fileext = ".fastq")
    recon_tmp <- tempfile(fileext = ".fastq")
    fraq_convert(src_path, src_tmp, nthreads = nthreads)
    fraq_convert(converted_path, recon_tmp, nthreads = nthreads)
    stopifnot(identical(readLines(src_tmp), readLines(recon_tmp)))
}

set.seed(1)
nthreads <- 1L

if (requireNamespace("ShortRead", quietly = TRUE)) {
    cat("Testing ShortRead import/export helpers...\n")

    # Single dataset round-trip ---------------------------------------------
    single_fastq <- tempfile(fileext = ".fastq")
    generate_random_fastq(
        single_fastq,
        n_reads = 24,
        read_length = 60,
        name_prefix = "srq_single_"
    )
    single_fraq <- tempfile(fileext = ".fraq")
    fraq_convert(single_fastq, single_fraq, nthreads = nthreads)
    single_srq <- fraq_export_shortreadq(single_fraq, nthreads = nthreads)
    stopifnot(inherits(single_srq, "ShortReadQ"))
    stopifnot(length(single_srq) == 24L)

    single_fastq_out <- tempfile(fileext = ".fastq.gz")
    fraq_import_shortreadq(single_srq, single_fastq_out, nthreads = nthreads)
    compare_fastq_records(single_fastq, single_fastq_out, nthreads = nthreads)

    # Paired dataset round-trip ---------------------------------------------
    paired_fastq <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq"))
    generate_random_fastq(
        paired_fastq[1],
        n_reads = 18,
        read_length = 55,
        name_prefix = "srq_pair1_"
    )
    generate_random_fastq(
        paired_fastq[2],
        n_reads = 18,
        read_length = 45,
        name_prefix = "srq_pair2_"
    )
    paired_mem <- c(tempfile(fileext = ".mem"), tempfile(fileext = ".mem"))
    fraq_convert(paired_fastq, paired_mem, nthreads = nthreads)
    paired_srq <- fraq_export_shortreadq(paired_mem, nthreads = nthreads)
    stopifnot(is.list(paired_srq))
    stopifnot(length(paired_srq) == 2L)
    stopifnot(all(vapply(paired_srq, inherits, logical(1), "ShortReadQ")))
    stopifnot(all(vapply(paired_srq, length, integer(1)) == 18L))

    paired_out <- c(tempfile(fileext = ".fastq"), tempfile(fileext = ".fastq.gz"))
    fraq_import_shortreadq(paired_srq, paired_out, nthreads = nthreads)
    compare_fastq_records(paired_fastq[1], paired_out[1], nthreads = nthreads)
    compare_fastq_records(paired_fastq[2], paired_out[2], nthreads = nthreads)

    cat("ShortRead helper tests completed successfully.\n")
} else {
    cat("Skipping ShortRead import/export helper tests; ShortRead is not installed.\n")
}

cat("Testing fraq_align Biostrings integration...\n")
align_queries <- Biostrings::DNAStringSet(c("ACGTNT", "TTT"))
align_targets <- Biostrings::DNAStringSet(c("ACGTAT", "GTTTGG"))
align_res <- fraq_align(
    align_queries,
    align_targets,
    max_distance = 2L,
    boundary = "contains",
    distance_metric = "lv"
)
stopifnot(inherits(align_res$query, "BStringSet"))
stopifnot(inherits(align_res$target, "BStringSet"))
stopifnot(identical(names(align_res)[1:2], c("query", "target")))
stopifnot(all(as.character(align_res$query) == as.character(align_queries)))
stopifnot(all(as.character(align_res$target) == as.character(align_targets)))
stopifnot(all(align_res$distance <= 2L))

align_hm <- fraq_align(
    "ACGT",
    "ACGG",
    max_distance = 1L,
    boundary = "global",
    distance_metric = "hm"
)
stopifnot(align_hm$distance[1] == 1L)
stopifnot(inherits(align_hm$query, "BStringSet"))
stopifnot(inherits(align_hm$target, "BStringSet"))
stopifnot(align_hm$start[1] == 0L)
stopifnot(align_hm$end[1] == 3L)

cat("Util function tests completed successfully.\n")
