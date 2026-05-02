# fraq 1.1.1

- Move ShortRead from Imports to Suggests and check for it at runtime in
  ShortRead bridge helpers.
- Add `fraq_run_r()` for block-oriented R kernels that run on the main R
  thread while file I/O and compression continue on `io_threads` TBB worker
  threads.
- Update bundled zstd to 1.5.7.
- Add TBB flow graph compatibility helpers for classic TBB and oneTBB APIs.

# fraq 0.99.2 (2026-02-12)
- Add user interrupt checks
- Add package man page
- Add installation instructions in vignette

# fraq 0.99.0 (2025-11-11)

- Initial Bioconductor submission.
