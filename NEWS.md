# fraq 1.1.1

- Move ShortRead from Imports to Suggests and check for it at runtime in
  ShortRead bridge helpers.
- Add `fraq_run_r()` for block-oriented R kernels that run on the main R
  thread while file I/O and compression can continue on `nthreads` TBB worker
  threads.
- Add serial fallback paths for C++ and R kernels when `nthreads = 1`, and
  force the serial path inside forked R processes to avoid using TBB after
  `fork()`.
- Update bundled zstd to 1.5.7.

# fraq 0.99.2 (2026-02-12)
- Add user interrupt checks
- Add package man page
- Add installation instructions in vignette

# fraq 0.99.0 (2025-11-11)

- Initial Bioconductor submission.
