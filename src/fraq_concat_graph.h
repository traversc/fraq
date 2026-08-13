#ifndef _FRAQ_CONCAT_GRAPH_H_
#define _FRAQ_CONCAT_GRAPH_H_

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <tbb/global_control.h>
#include "fraq_defines.h"
#include "io/fraq_format.h"
#include "io/fraq_readers.h"
#include "io/tbb_flow_compat.h"
#include "io/fraq_writers.h"

namespace fraq_internal {
using fraq::Read;
using fraq::ends_with;
using tbb::flow::continue_msg;
using tbb::flow::function_node;
using tbb::flow::graph;
using tbb::flow::limiter_node;
using tbb::flow::make_edge;
using tbb::flow::serial;
using tbb::flow::sequencer_node;
using tbb::flow::unlimited;

struct ConcatReadsBlock {
  size_t index{0};
  std::vector<Read> reads;
};
using ConcatReadsBlockPtr = std::shared_ptr<ConcatReadsBlock>;

struct ConcatCompressedBlock {
  size_t index{0};
  CompressedReadBlock block;
};
using ConcatCompressedBlockPtr = std::shared_ptr<ConcatCompressedBlock>;

inline MultiReader make_reader(const std::string &path,
                               std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store) {
  if (ends_with(path, ".mem")) {
    return MultiReader(path, mem_store);
  }
  return MultiReader(path);
}

inline int pick_fastq_compress_level(const std::string &path, const fraq::FraqRunConfig &config) {
  if (ends_with(path, ".zst")) return config.zstd_compress_level;
  return config.gzip_compress_level; // plain writer ignores level
}

struct FastqConcatGraph {
  static constexpr size_t MAX_INFLIGHT_INDICES = 4;

  const fraq::FraqRunConfig config;
  MultiWriter writer;
  std::vector<MultiReader> readers;
  graph flow_graph;
  limiter_node<continue_msg> read_limiter;
  function_node<continue_msg, ConcatReadsBlockPtr> read_node;
  function_node<ConcatReadsBlockPtr, continue_msg> write_node;
  size_t current_reader{0};
  size_t next_index{0};
  bool exhausted{false};

  FastqConcatGraph(const std::vector<std::string> &input_files,
                       const fraq::FraqRunConfig &cfg,
                       std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store,
                       const std::string &output_file)
      : config(cfg),
        writer(output_file, pick_fastq_compress_level(output_file, cfg)),
        readers(),
        flow_graph(),
        read_limiter(flow_graph, MAX_INFLIGHT_INDICES),
        read_node(flow_graph, serial, [this](const continue_msg &) { return this->read_step(); }),
        write_node(flow_graph, serial, [this](const ConcatReadsBlockPtr &ptr) { return this->write_step(ptr); }) {
    readers.reserve(input_files.size());
    for (const auto &file : input_files) {
      readers.emplace_back(make_reader(file, mem_store));
    }
    make_edge(read_limiter, read_node);
    make_edge(read_node, write_node);
  }

  ConcatReadsBlockPtr read_step() {
    if (exhausted) {
      return {};
    }
    while (current_reader < readers.size()) {
      auto block_ptr = std::make_shared<ConcatReadsBlock>();
      block_ptr->index = next_index;
      block_ptr->reads.reserve(config.blocksize);
      readers[current_reader].read_block(block_ptr->reads, config.blocksize);
      if (block_ptr->reads.empty()) {
        ++current_reader;
        continue;
      }
      ++next_index;
      return block_ptr;
    }
    exhausted = true;
    return {};
  }

  continue_msg write_step(const ConcatReadsBlockPtr &ptr) {
    tbb_compat::decrementer(read_limiter).try_put(continue_msg{});
    if (ptr && !ptr->reads.empty()) {
      writer.write_block(ptr->reads);
      read_limiter.try_put(continue_msg{});
    }
    return continue_msg{};
  }

  void start() {
    for (size_t i = 0; i < MAX_INFLIGHT_INDICES; ++i) {
      read_limiter.try_put(continue_msg{});
    }
  }

  void wait() {
    flow_graph.wait_for_all();
  }

  void flush() {
    // FIFO writers buffer in memory instead of blocking when the pipe is full,
    // so bytes can still be pending here. close() would discard them, so drain
    // first, the same way FraqRunGraph::wait_and_flush() does.
    while (writer.has_pending()) {
      writer.flush();
    }
    writer.close();
  }
};

struct FraqConcatGraph {
  static constexpr size_t MAX_INFLIGHT_INDICES = 4;

  const fraq::FraqRunConfig config;
  std::shared_ptr<FraqfMultiWriter> writer;
  std::vector<MultiReader> readers;
  graph flow_graph;
  limiter_node<continue_msg> read_limiter;
  function_node<continue_msg, continue_msg> read_node;
  function_node<ConcatReadsBlockPtr, ConcatCompressedBlockPtr> compress_node;
  sequencer_node<ConcatCompressedBlockPtr> order_node;
  function_node<ConcatCompressedBlockPtr, continue_msg> write_node;
  std::vector<Read> scratch;
  size_t current_reader{0};
  size_t next_index{0};
  bool exhausted{false};

  FraqConcatGraph(const std::vector<std::string> &input_files,
                      const fraq::FraqRunConfig &cfg,
                      std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store,
                      const std::string &output_file)
      : config(cfg),
        writer(ends_with(output_file, ".fraq")
                   ? std::make_shared<FraqfMultiWriter>(output_file)
                   : std::make_shared<FraqfMultiWriter>(output_file, mem_store)),
        readers(),
        flow_graph(),
        read_limiter(flow_graph, MAX_INFLIGHT_INDICES),
        read_node(flow_graph, serial, [this](const continue_msg &) { return this->read_step(); }),
        compress_node(flow_graph,
                      unlimited,
                      [this](const ConcatReadsBlockPtr &ptr) { return this->compress_step(ptr); }),
        order_node(flow_graph,
                   [](const ConcatCompressedBlockPtr &ptr) { return ptr->index; }),
        write_node(flow_graph,
                   serial,
                   [this](ConcatCompressedBlockPtr ptr) { return this->write_step(std::move(ptr)); }),
        scratch() {
    readers.reserve(input_files.size());
    for (const auto &file : input_files) {
      readers.emplace_back(make_reader(file, mem_store));
    }
    scratch.reserve(config.blocksize);
    make_edge(read_limiter, read_node);
    make_edge(compress_node, order_node);
    make_edge(order_node, write_node);
  }

  continue_msg read_step() {
    if (exhausted) {
      tbb_compat::decrementer(read_limiter).try_put(continue_msg{});
      return continue_msg{};
    }
    while (current_reader < readers.size()) {
      auto block_ptr = std::make_shared<ConcatReadsBlock>();
      block_ptr->index = next_index;
      block_ptr->reads.reserve(config.blocksize);
      while (block_ptr->reads.size() < config.blocksize && current_reader < readers.size()) {
        const size_t remaining = config.blocksize - block_ptr->reads.size();
        readers[current_reader].read_block(scratch, remaining);
        if (scratch.empty()) {
          ++current_reader;
          continue;
        }
        const bool reader_drained = scratch.size() < remaining;
        block_ptr->reads.insert(block_ptr->reads.end(),
                                std::make_move_iterator(scratch.begin()),
                                std::make_move_iterator(scratch.end()));
        scratch.clear();
        if (reader_drained) {
          ++current_reader;
        }
      }
      if (block_ptr->reads.empty()) {
        continue;
      }
      ++next_index;
      compress_node.try_put(block_ptr);
      return continue_msg{};
    }
    exhausted = true;
    tbb_compat::decrementer(read_limiter).try_put(continue_msg{});
    return continue_msg{};
  }

  ConcatCompressedBlockPtr compress_step(const ConcatReadsBlockPtr &ptr) {
    if (!ptr) return {};
    auto compressed_ptr = std::make_shared<ConcatCompressedBlock>();
    compressed_ptr->index = ptr->index;
    compressed_ptr->block = compress_reads(ptr->reads, config.fraq_compress_level);
    return compressed_ptr;
  }

  continue_msg write_step(ConcatCompressedBlockPtr ptr) {
    tbb_compat::decrementer(read_limiter).try_put(continue_msg{});
    if (ptr) {
      writer->write_compressed_block(std::move(ptr->block));
      read_limiter.try_put(continue_msg{});
    }
    return continue_msg{};
  }

  void start() {
    for (size_t i = 0; i < MAX_INFLIGHT_INDICES; ++i) {
      read_limiter.try_put(continue_msg{});
    }
  }

  void wait() {
    flow_graph.wait_for_all();
  }

  void flush() {
    writer->close();
  }
};

} // namespace fraq_internal

#endif // _FRAQ_CONCAT_GRAPH_H_
