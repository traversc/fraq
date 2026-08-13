# fraq 1.0.2

- Fix `fraq_concat()` discarding output written to a `.fifo` path when the
  reader is not keeping up. Buffered bytes are now drained before the writer
  is closed, matching what `fraq_run()` already did.

# fraq 1.0.1

- Add TBB flow graph compatibility helpers so the package builds against
  both classic TBB and oneTBB (fixes build failure caused by an
  RcppParallel/TBB update that removed `tbb::flow::source_node` and made
  `limiter_node::decrement` private).
- Fix `.zst` outputs gaining a spurious trailing empty frame: `ZstdWriter`
  now flushes with `ZSTD_flushStream` and only ends the frame in `close()`.
- Close `.fraq`/`.mem` writers explicitly on the success path, so write
  errors at close time are reported instead of silently discarded.

# fraq 0.99.2 (2026-02-12)
- Add user interrupt checks
- Add package man page
- Add installation instructions in vignette

# fraq 0.99.0 (2025-11-11)

- Initial Bioconductor submission.
