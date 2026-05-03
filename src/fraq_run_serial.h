#ifndef _FRAQ_RUN_SERIAL_H_
#define _FRAQ_RUN_SERIAL_H_

#include <algorithm>
#include <cstdio>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "fraq_run_graph_r.h"

namespace fraq_internal {

inline MultiReader make_serial_reader(
    const std::string &path,
    std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store) {
  if (fraq::ends_with(path, ".mem")) {
    return MultiReader(path, mem_store);
  }
  return MultiReader(path);
}

struct SerialFraqSink {
  std::unique_ptr<FraqfMultiWriter> writer;
  std::vector<Read> current_block;

  SerialFraqSink(const std::string &path,
                 const size_t blocksize,
                 std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store)
      : writer(fraq::ends_with(path, ".mem")
                   ? std::make_unique<FraqfMultiWriter>(path, mem_store)
                   : std::make_unique<FraqfMultiWriter>(path)),
        current_block() {
    current_block.reserve(blocksize);
  }
};

struct SerialOutputRouter {
  const FraqRunConfig &config;
  std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store;
  std::unordered_map<std::string, std::unique_ptr<MultiWriter>> fastq_writers;
  std::unordered_map<std::string, std::unique_ptr<SerialFraqSink>> fraq_sinks;
  bool closed{false};

  SerialOutputRouter(
      const FraqRunConfig &config,
      std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store)
      : config(config),
        mem_store(mem_store),
        fastq_writers(),
        fraq_sinks() {}

  ~SerialOutputRouter() {
    try {
      close_all();
    } catch (...) {
    }
  }

  void write_processed(const ProcessedBlockPtr &p) {
    if (!p) {
      return;
    }
    for (auto &kv : p->bins) {
      const std::string &path = kv.first;
      std::vector<Read> &batch = kv.second;
      if (batch.empty()) {
        continue;
      }
      if (fraq::ends_with(path, ".fraq") || fraq::ends_with(path, ".mem")) {
        write_fraq_batch(path, batch);
      } else {
        MultiWriter &writer = get_fastq_writer(path);
        writer.write_block(batch);
      }
    }
  }

  void finish() {
    for (auto &kv : fraq_sinks) {
      flush_fraq_sink(*kv.second);
    }

    bool pending = true;
    while (pending) {
      pending = false;
      for (auto &kv : fastq_writers) {
        MultiWriter &writer = *kv.second;
        if (writer.has_pending()) {
          writer.flush();
          if (writer.has_pending()) {
            pending = true;
          }
        }
      }
    }
    close_all();
  }

  void cleanup_outputs() {
    if (closed) {
      return;
    }
    for (auto &kv : fastq_writers) {
      const std::string &path = kv.first;
      try {
        kv.second->close();
      } catch (...) {
      }
      if (!fraq::ends_with(path, ".fifo")) {
        std::remove(path.c_str());
      }
    }
    for (auto &kv : fraq_sinks) {
      const std::string &path = kv.first;
      try {
        kv.second->writer->close();
      } catch (...) {
      }
      if (fraq::ends_with(path, ".mem")) {
        mem_store.erase(path);
      } else {
        std::remove(path.c_str());
      }
    }
    closed = true;
  }

private:
  MultiWriter &get_fastq_writer(const std::string &path) {
    auto it = fastq_writers.find(path);
    if (it != fastq_writers.end()) {
      return *it->second;
    }
    const int compress_level = fraq::ends_with(path, ".gz") ? config.gzip_compress_level :
                               fraq::ends_with(path, ".zst") ? config.zstd_compress_level :
                               1;
    auto writer = std::make_unique<MultiWriter>(path, compress_level);
    MultiWriter *ptr = writer.get();
    fastq_writers.emplace(path, std::move(writer));
    return *ptr;
  }

  SerialFraqSink &get_fraq_sink(const std::string &path) {
    auto it = fraq_sinks.find(path);
    if (it != fraq_sinks.end()) {
      return *it->second;
    }
    auto sink = std::make_unique<SerialFraqSink>(path, config.blocksize, mem_store);
    SerialFraqSink *ptr = sink.get();
    fraq_sinks.emplace(path, std::move(sink));
    return *ptr;
  }

  void write_fraq_batch(const std::string &path, std::vector<Read> &batch) {
    SerialFraqSink &sink = get_fraq_sink(path);
    for (auto &read : batch) {
      sink.current_block.push_back(std::move(read));
      if (sink.current_block.size() >= config.blocksize) {
        flush_fraq_sink(sink);
      }
    }
  }

  void flush_fraq_sink(SerialFraqSink &sink) {
    if (sink.current_block.empty()) {
      return;
    }
    CompressedReadBlock compressed = compress_reads(sink.current_block, config.fraq_compress_level);
    sink.writer->write_compressed_block(std::move(compressed));
    sink.current_block.clear();
  }

  void close_all() {
    if (closed) {
      return;
    }
    for (auto &kv : fastq_writers) {
      kv.second->close();
    }
    for (auto &kv : fraq_sinks) {
      kv.second->writer->close();
    }
    closed = true;
  }
};

inline std::vector<size_t> serial_excess_readers(const std::vector<size_t> &read_counts) {
  std::vector<size_t> excess_indices;
  if (read_counts.empty()) {
    return excess_indices;
  }
  const size_t min_reads = *std::min_element(read_counts.begin(), read_counts.end());
  for (size_t i = 0; i < read_counts.size(); ++i) {
    if (read_counts[i] > min_reads) {
      excess_indices.push_back(i);
    }
  }
  return excess_indices;
}

inline bool serial_interrupt_requested(const FraqRunConfig &config) {
  return config.interrupt_fn && config.interrupt_fn(config.interrupt_ctx);
}

inline bool serial_read_joined_block(std::vector<MultiReader> &readers,
                                     const size_t block_index,
                                     const FraqRunConfig &config,
                                     std::vector<size_t> &read_counts,
                                     BlockPtrVec &joined) {
  joined = BlockPtrVec(readers.size());
  bool any_reads = false;
  for (size_t ri = 0; ri < readers.size(); ++ri) {
    auto block = std::make_shared<Block>();
    block->index = block_index;
    block->reader_index = ri;
    readers[ri].read_block(block->reads, config.blocksize);
    read_counts[ri] += block->reads.size();
    if (!block->reads.empty()) {
      joined.insert(block, ri);
      any_reads = true;
    }
  }
  return any_reads;
}

inline std::vector<size_t> serial_fraq_run(
    const std::vector<std::string> &input_files,
    fraq::process_task_t process_kernel,
    const FraqRunConfig &config,
    std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store) {
  if (input_files.empty()) {
    throw std::runtime_error("fraq_run requires at least one input file");
  }

  std::vector<MultiReader> readers;
  readers.reserve(input_files.size());
  for (const auto &file : input_files) {
    readers.emplace_back(make_serial_reader(file, mem_store));
  }

  std::vector<size_t> read_counts(input_files.size(), 0);
  SerialOutputRouter router(config, mem_store);

  try {
    for (size_t block_index = 0;; ++block_index) {
      if (serial_interrupt_requested(config)) {
        router.cleanup_outputs();
        return serial_excess_readers(read_counts);
      }
      const size_t first_index = block_index * config.blocksize;
      if (config.limit > 0 && first_index >= config.limit) {
        break;
      }

      BlockPtrVec joined;
      const bool any_reads = serial_read_joined_block(readers, block_index, config, read_counts, joined);
      if (!any_reads || !joined.full()) {
        break;
      }

      auto processed = std::make_shared<ProcessedBlock>();
      processed->index = block_index;
      const size_t nreads = joined.min_reads();
      for (size_t i = 0; i < nreads; ++i) {
        const size_t id = first_index + i;
        if (config.limit > 0 && id >= config.limit) {
          break;
        }
        std::vector<Read> input_reads = joined.move_reads(i);
        auto routed = process_kernel(std::move(input_reads), id);
        for (auto &pr : routed) {
          if (pr.first == fraq::NO_OUTPUT) {
            continue;
          }
          processed->bins[pr.first].push_back(std::move(pr.second));
        }
      }
      router.write_processed(processed);
    }
    router.finish();
  } catch (...) {
    router.cleanup_outputs();
    throw;
  }

  return serial_excess_readers(read_counts);
}

inline std::vector<size_t> serial_fraq_run_r(
    const std::vector<std::string> &input_files,
    Rcpp::Function &kernel,
    const FraqRunConfig &config,
    std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store) {
  if (input_files.empty()) {
    throw std::runtime_error("fraq_run_r requires at least one input file");
  }

  std::vector<MultiReader> readers;
  readers.reserve(input_files.size());
  for (const auto &file : input_files) {
    readers.emplace_back(make_serial_reader(file, mem_store));
  }

  std::vector<size_t> read_counts(input_files.size(), 0);
  SerialOutputRouter router(config, mem_store);

  try {
    for (size_t block_index = 0;; ++block_index) {
      if (serial_interrupt_requested(config)) {
        router.cleanup_outputs();
        return serial_excess_readers(read_counts);
      }
      Rcpp::checkUserInterrupt();
      const size_t first_index = block_index * config.blocksize;
      if (config.limit > 0 && first_index >= config.limit) {
        break;
      }

      BlockPtrVec joined;
      const bool any_reads = serial_read_joined_block(readers, block_index, config, read_counts, joined);
      if (!any_reads || !joined.full()) {
        break;
      }

      ProcessedBlockPtr processed = r_kernel_process_block(joined, kernel, config);
      router.write_processed(processed);
    }
    router.finish();
  } catch (...) {
    router.cleanup_outputs();
    throw;
  }

  return serial_excess_readers(read_counts);
}

} // namespace fraq_internal

#endif // _FRAQ_RUN_SERIAL_H_
