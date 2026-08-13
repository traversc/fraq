# fraq 1.1.3

- Fix `fraq_concat()` discarding output written to a `.fifo` path when the
  reader is not keeping up. Buffered bytes are now drained before the writer
  is closed, matching what `fraq_run()` already did.
- Fix a file descriptor leak when closing a `.fifo` writer failed.
- Detect forked R processes with `RcppParallel::isProcessForkedChild()`
  instead of probing `parallel::isChild()`, so the serial fallback triggers
  regardless of which package performed the fork. This raises the minimum
  RcppParallel version to 6.1.1.
- `PlainWriter` and `FraqfFileWriter` now tolerate being closed twice, like
  the gzip and zstd writers already did.

# fraq 1.1.2

- Fix `.zst` outputs gaining a spurious trailing empty frame: `ZstdWriter`
  now flushes with `ZSTD_flushStream` and only ends the frame in `close()`.
- Close `.fraq`/`.mem` writers explicitly on the success path, so write
  errors at close time are reported instead of silently discarded.

# fraq 1.1.1

- Move ShortRead from Imports to Suggests and check for it at runtime in
  ShortRead bridge helpers.
- Add TBB flow graph compatibility helpers so the package builds against
  both classic TBB and oneTBB (fixes a build failure caused by an
  RcppParallel/TBB update that removed `tbb::flow::source_node` and made
  `limiter_node::decrement` private).
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
