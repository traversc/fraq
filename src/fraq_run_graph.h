#ifndef _FRAQ_RUN_GRAPH_H
#define _FRAQ_RUN_GRAPH_H

#include <zlib.h>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include <utility>
#include <limits>
#include <atomic>
#include <stdexcept>
#include <sstream>
#include <memory>
#include <tbb/flow_graph.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/task_group.h>
#include "io/fraq_readers.h"
#include "io/fraq_writers.h"
#include "fraq_defines.h"

namespace fraq_internal {
using namespace fraq;

using tbb::flow::graph;
using tbb::flow::source_node;
using tbb::flow::function_node;
using tbb::flow::sequencer_node;
using tbb::flow::continue_msg;
using tbb::flow::serial;
using tbb::flow::unlimited;
using tbb::flow::limiter_node;

// atomic min helper function
template<typename T>
inline void atomic_min(std::atomic<T>& a, T v) {
    T old = a.load(std::memory_order_relaxed);
    while (old > v && !a.compare_exchange_weak(old, v, std::memory_order_relaxed)) {}
}

struct Block {
  std::vector<Read> reads;
  size_t index; // block index
  size_t reader_index;
};
using BlockPtr = std::shared_ptr<Block>;

// blocks from different fastqs with the same index
// e.g. in paired end sequencing blocks is length two with the same index, one for each mate
struct BlockPtrVec {
  std::vector<BlockPtr> blocks;
  size_t added_blocks;
  BlockPtrVec() : added_blocks(0) {}
  BlockPtrVec(const size_t n) : blocks(n), added_blocks(0) {}
  void insert(const BlockPtr &b, const size_t index) {
    auto& slot = blocks.at(index);
    if(slot) throw std::runtime_error("Internal error: duplicate Block for this reader/index");
    slot = b;
    ++added_blocks;
  }
  size_t size() const {
    return blocks.size();
  }
  bool full() const {
    return added_blocks == blocks.size();
  }
  void resize(const size_t n) {
    blocks.resize(n);
  }
  size_t min_reads() const {
    size_t min_size = std::numeric_limits<size_t>::max();
    for (const auto &b : blocks) {
      if (b && b->reads.size() < min_size) {
        min_size = b->reads.size();
      }
    }
    return min_size;
  }
  // destructively move reads to output
  std::vector<Read> move_reads(const size_t i) {
    std::vector<Read> output_reads;
    output_reads.reserve(blocks.size());
    for (auto &b : blocks) {
      output_reads.push_back(std::move(b->reads[i]));
    }
    return output_reads;
  }
};

// Processed block with target filenames for each mate
struct ProcessedBlock {
  std::unordered_map<std::string, std::vector<Read>> bins; // fname -> reads
  size_t index;
};
using ProcessedBlockPtr = std::shared_ptr<ProcessedBlock>;

// Writer state for demultiplexing
struct DemuxWriter {
  std::vector<Read> current_block; // accessed only single threaded by demux_node
  tbb::concurrent_queue<std::vector<Read>> queue; // may be accessed concurrently
  // in_flight is a semaphore to avoid over-scheduling writes to file
  // Note: we cannot use std::mutex as mutex must be locked/unlocked from the same thread, 
  // whereas the scheduler (demux_node) can be a different thread than the writer_node where it would be unlocked
  std::atomic<bool> in_flight{false};
  MultiWriter writer;
  DemuxWriter(const std::string &fname, const size_t blocksize, const int compress_level) : writer(fname, compress_level) {
    current_block.reserve(blocksize);
  }
  size_t current_block_size() const {
    return current_block.size();
  }
  size_t add_read(Read && read) { // destructive add
    current_block.push_back(std::move(read));
    return current_block.size();
  }
  void push_current_block() {
    if (current_block.empty()) return; // nothing to write
    std::vector<Read> new_current_block;
    new_current_block.reserve(current_block.size());
    queue.push(std::move(current_block));
    current_block = std::move(new_current_block); // clear current_block and replace it with an empty vector
  }
};
using DemuxWriterPtr = std::shared_ptr<DemuxWriter>;

// fraq format (fraqf) demux writer
// This differs from FastqDemuxWriter in that block compression is multithreaded so is more complex
// FraqfDemuxWriter contains all writer info and compressed blocks
// But we need to release read blocks so that they can be compressed async while work is being done on an instance
struct FraqfDemuxWriter;
struct FraqfOutputBlock {
  std::vector<Read> reads;
  size_t index;
  FraqfDemuxWriter * w;
};
using FraqfOutputBlockPtr = std::shared_ptr<FraqfOutputBlock>;
struct FraqfDemuxWriter {
  // accessed only single threaded by demux_node
  std::vector<Read> current_block;
  size_t current_block_index{0};
  tbb::concurrent_unordered_map<size_t, CompressedReadBlock> compressed_blocks;

  // multithreaded access by demux_node and fraqf_compressor_node and fraqf_writer_node
  // in_flight is a semaphore to avoid over-scheduling writes to file
  std::atomic<bool> in_flight{false};
  size_t next_expected_index{0}; // next expected index to write
  FraqfMultiWriter writer;
  FraqfDemuxWriter(const std::string &fname, const size_t blocksize) : writer(fname) {
    current_block.reserve(blocksize);
  } // .fraq
  FraqfDemuxWriter(const std::string &mem_key, const size_t blocksize, std::unordered_map<std::string, std::vector<CompressedReadBlock>> & store) : writer(mem_key, store) {
    current_block.reserve(blocksize);
  } // .mem
  size_t current_block_size() const {
    return current_block.size();
  }
  size_t add_read(Read && read) { // destructive add
    current_block.push_back(std::move(read));
    return current_block.size();
  }
  FraqfOutputBlockPtr release_current_block() {
    std::vector<Read> new_current_block;
    new_current_block.reserve(current_block.size());
    auto block = std::make_shared<FraqfOutputBlock>();
    block->reads = std::move(current_block);
    block->index = current_block_index;
    block->w = this;
    current_block = std::move(new_current_block);
    ++current_block_index;
    return block;
  }
  // r-value reference prevents copying of compressed_block
  void add_compressed_block(size_t index, CompressedReadBlock &&compressed_block) {
    compressed_blocks.emplace(index, std::move(compressed_block));
  }
  void write_next_compressed_blocks() {
    while (write_next_compressed_block()) {}
  }
  private:
  bool write_next_compressed_block() {
    auto it = compressed_blocks.find(next_expected_index);
    if (it == compressed_blocks.end()) {
      return false;
    }
    writer.write_compressed_block(std::move(it->second));
    compressed_blocks.unsafe_erase(it);
    ++next_expected_index;
    return true;
  }
};
using FraqfDemuxWriterPtr = std::shared_ptr<FraqfDemuxWriter>;

struct FraqRunGraph {

  // Max number of concurrently in-flight joined indices
  // Controls backpressure on the primary reader
  static constexpr size_t MAX_INFLIGHT_INDICES = 4;

  graph flow_graph;
  const FraqRunConfig config;
  std::unordered_map<std::string, std::vector<CompressedReadBlock>> & mem_store; // external memory store for .mem files
  
  // Reader nodes
  // block_index_limit will be atomicly updated if one of the readers reaches EOF
  // Subsequently, all readers will stop reading when their block index reaches this limit
  std::vector<MultiReader> readers;
  std::vector<size_t> current_index;
  std::vector<size_t> read_count_per_reader;
  std::atomic<size_t> block_index_limit{std::numeric_limits<size_t>::max()}; // last valid block index (inclusive)
  // Primary and secondary reader nodes (lockstep scheduling)
  limiter_node<continue_msg> primary_limiter;           // gates how many indices can be in flight
  function_node<continue_msg, BlockPtr> primary_reader_node; // drives reader_index = 0
  std::vector< function_node<continue_msg, BlockPtr> > secondary_reader_nodes; // drives reader_index = 1..N-1

  // Mate joining node, joins read blocks keyed on index, stateful and serial
  // Use unordered_map to bound memory to in-flight indices only
  std::unordered_map<size_t, BlockPtrVec> read_blocks;
  function_node<BlockPtr, continue_msg> joiner_node;

  // Kernel processing node
  process_task_t process_kernel; // kernel
  function_node<BlockPtrVec, ProcessedBlockPtr> process_node;

  // Reorder node, reorders processed blocks by index and passes it to demux node
  sequencer_node<ProcessedBlockPtr> reorder_node;

  // Demux node adds reads to the appropriate queue (DemuxWriter) based on output filenames
  // Constructs a DemuxWriter on the fly if not already present
  // sends a job signal to writer_node when a DemuxWriter's queue reaches BLOCKSIZE
  std::unordered_map<std::string, DemuxWriterPtr> fastq_writer_map;
  std::unordered_map<std::string, FraqfDemuxWriterPtr> fraqf_writer_map;
  function_node<ProcessedBlockPtr, continue_msg> demux_node;
  // Global write queue for FASTQ writers to enable explicit draining
  tbb::concurrent_queue<DemuxWriter*> fastq_write_queue;
  
  // Writer for fastq format (.fq, .gz, .zst)
  function_node<continue_msg, continue_msg> fastq_writer_node;

  // Writer node for fraqf format (.fraq, .mem)
  // Compression is multithreaded at the block level, so needs a separate node
  tbb::concurrent_queue<FraqfOutputBlockPtr> fraqf_compress_queue;
  tbb::concurrent_queue<FraqfDemuxWriter*> fraqf_write_queue;
  function_node<continue_msg, continue_msg> fraqf_compressor_node;
  function_node<continue_msg, continue_msg> fraqf_writer_node;

  FraqRunGraph(const std::vector<std::string> &input_files,
                 process_task_t process_kernel,
                 const FraqRunConfig config,
                 std::unordered_map<std::string, std::vector<CompressedReadBlock>> & mem_store) :
    config(config),
    mem_store(mem_store),
    readers(readers_init(input_files)),
    current_index(input_files.size(), 0),
    read_count_per_reader(input_files.size(), 0),
    primary_limiter(flow_graph, MAX_INFLIGHT_INDICES),
    primary_reader_node(
      flow_graph,
      serial,
      [this](const continue_msg&) -> BlockPtr { return this->primary_reader_step(); }
    ),
    joiner_node(flow_graph, serial,
      [this](const BlockPtr &b) { return this->joiner_node_body(b); }
    ),
    process_kernel(process_kernel),
    process_node(flow_graph, config.serial_kernel ? serial : unlimited, 
      [this](BlockPtrVec b) { return this->process_node_body(b); }
    ),
    reorder_node(flow_graph,
      [](const ProcessedBlockPtr &p) { return p->index; }
    ),
    demux_node(flow_graph, serial,
      [this](const ProcessedBlockPtr &p) { return this->demux_node_body(p); }
    ),
    fastq_writer_node(flow_graph, unlimited,
      [this](const continue_msg&) { return this->fastq_writer_node_body(); }
    ),
    fraqf_compressor_node(flow_graph, unlimited,
      [this](const continue_msg&) { return this->fraqf_compressor_node_body(); }
    ),
    fraqf_writer_node(flow_graph, unlimited,
      [this](const continue_msg&) { return this->fraqf_writer_node_body(); }
    )
    {
    // Initialize secondary reader nodes for reader_index = 1..N-1
    if (readers.size() > 1) {
      secondary_reader_nodes.reserve(readers.size() - 1);
      for (size_t ri = 1; ri < readers.size(); ++ri) {
        secondary_reader_nodes.emplace_back(
          flow_graph,
          serial,
          [this, ri](const continue_msg&) -> BlockPtr {
            return this->reader_node_body(ri);
          }
        );
      }
    }

    // Wire limiter -> primary reader, and readers to the joiner
    make_edge(primary_limiter, primary_reader_node);
    make_edge(primary_reader_node, joiner_node);
    for (auto &fn : secondary_reader_nodes) {
      make_edge(fn, joiner_node);
    }
    make_edge(process_node, reorder_node);
    make_edge(reorder_node, demux_node);
    // Seed limiter with initial K requests
    for(size_t k = 0; k < MAX_INFLIGHT_INDICES; ++k) {
      primary_limiter.try_put(continue_msg());
    }
  }
  // Debug summary removed.

  // final synchronization and making sure all writers and buffers are flushed
  void wait_and_flush() {
    // wait for all data processing to finish
    flow_graph.wait_for_all();

    // Concurrent flush of FRAQ writers
    for (auto& kv : fraqf_writer_map) {
      FraqfDemuxWriter* w = kv.second.get();
      if (w->current_block_size() > 0) {
        FraqfOutputBlockPtr ob = w->release_current_block();
        fraqf_compress_queue.push(ob);
        fraqf_compressor_node.try_put(continue_msg{});
      }
    }
    // Ensure FRAQ finishes by calling in to compressor--(worksteal)-->fraq_writer, single-threaded
    this->fraqf_compressor_node_body();

    // Concurrent flush of fastq writers
    for (auto& kv : fastq_writer_map) {
      DemuxWriter* w = kv.second.get();
      if (w->current_block_size() > 0) {
        w->push_current_block();
      }
      if (!w->in_flight.exchange(true, std::memory_order_acq_rel)) {
        fastq_write_queue.push(w);
        fastq_writer_node.try_put(continue_msg{});
      }
    }
    // finish all concurrent flushes
    flow_graph.wait_for_all();

    // fifo writers may still have buffered writes, because we didn't want to block when a fifo was full (which may cause deadlock)
    // Since we are at the end, we will clear all fifo writers and block by calling within a while loop
    // Its important to do this in a round robin fashion, in case one an external program is blocking one fifo writer because it expects output from another
    bool pending = true;
    while (pending) {
      pending = false;
      for (auto& kv : fastq_writer_map) {
        DemuxWriter* w = kv.second.get();
        if (w->writer.has_pending()) {
          w->writer.flush();
          if (w->writer.has_pending()) {
            pending = true;
          }
        }
      }
    }
    for (auto& kv : fastq_writer_map) {
      kv.second->writer.close();
    }
  }

  std::vector<size_t> check_reader_balance() const {
    std::vector<size_t> excess_indices;
    if (read_count_per_reader.empty()) return excess_indices;
    size_t min_reads = *std::min_element(read_count_per_reader.begin(), read_count_per_reader.end());
    for (size_t i = 0; i < read_count_per_reader.size(); ++i) {
      if (read_count_per_reader[i] > min_reads) {
        excess_indices.push_back(i);
      }
    }
    return excess_indices;
  }

  
  // initialization helpers
  std::vector<MultiReader> readers_init(const std::vector<std::string> &input_files) {
    std::vector<MultiReader> readers;
    readers.reserve(input_files.size());
    for (const auto &file : input_files) {
      if(fraq::ends_with(file, ".mem")) {
        readers.emplace_back(file, mem_store);
      } else {
        readers.emplace_back(file);
      }
    }
    return readers;
  }
  // node bodies
  // Primary reader (reader_index = 0). Reads one block and triggers secondary readers.
  BlockPtr primary_reader_step() {
      // Schedule secondary readers for this step BEFORE reading primary,
      // so other threads can begin work immediately.
      for (auto &fn : secondary_reader_nodes) {
        fn.try_put(continue_msg{});
      }
      return reader_node_body(0);
  }

  // ri = reader index, 0 for primary, 1..N-1 for secondary
  BlockPtr reader_node_body(const size_t ri) {
      size_t& idx = current_index[ri];
      const size_t limit = block_index_limit.load(std::memory_order_acquire);
      auto b = std::make_shared<Block>();
      b->index = idx;
      b->reader_index = ri;
      if (idx > limit) {
          // Past limit: emit empty block which joiner ignores
          b->reads.clear();
          return b;
      }
      readers[ri].read_block(b->reads, config.blocksize);
      read_count_per_reader[ri] += b->reads.size();
      if (b->reads.size() < config.blocksize) {
          atomic_min(block_index_limit, idx);
      }
      ++idx;
      return b;
  }

  continue_msg joiner_node_body(const BlockPtr &b) {
    // Ignore empty blocks (used to drain pending limiter tokens past EOF)
    if (!b || b->reads.empty()) {
      return continue_msg();
    }
    size_t idx = b->index;
    // Add block to the read_blocks map, keyed by index
    auto it = read_blocks.find(idx);
    if (it == read_blocks.end()) {
      BlockPtrVec bpv;
      bpv.resize(readers.size());
      auto ins = read_blocks.emplace(idx, std::move(bpv));
      it = ins.first;
      
    }
    it->second.insert(b, b->reader_index);
    
    // If all readers have added their blocks, process the block
    if(it->second.full()) {
      BlockPtrVec to_process = std::move(it->second);
      read_blocks.erase(it);
      process_node.try_put(std::move(to_process));
    }
    return continue_msg();
  }
  
  ProcessedBlockPtr process_node_body(BlockPtrVec b) {
    auto p = std::make_shared<ProcessedBlock>();
    const size_t nreads = b.min_reads();
    p->index = b.blocks[0]->index;

    for (size_t i = 0; i < nreads; ++i) {
      std::vector<Read> input_reads = b.move_reads(i);
      const size_t id = p->index * config.blocksize + i;

      // kernel returns vector<pair<fname, Read>>
      auto routed = process_kernel(input_reads, id);
      for (auto& pr : routed) {
        if (pr.first == NO_OUTPUT) continue;
        p->bins[pr.first].push_back(std::move(pr.second));
      }
    }
    
    return p;
  }

  continue_msg demux_node_body(const ProcessedBlockPtr& p) {
    
    for (auto& kv : p->bins) {
      const std::string& fname = kv.first;
      std::vector<Read>& batch  = kv.second;

      if (ends_with(fname, ".fraq") || ends_with(fname, ".mem")) {
        FraqfDemuxWriterPtr w;
        auto itf = fraqf_writer_map.find(fname);
        if (itf != fraqf_writer_map.end()) {
          w = itf->second;
        } else {
          w = ends_with(fname, ".mem")
            ? std::make_shared<FraqfDemuxWriter>(fname, config.blocksize, mem_store)
            : std::make_shared<FraqfDemuxWriter>(fname, config.blocksize);
          fraqf_writer_map[fname] = w;
        }
        for (auto& r : batch) {
          size_t sz = w->add_read(std::move(r));
          if (sz >= config.blocksize) {
            auto ob = w->release_current_block();
            if (ob) {
              fraqf_compress_queue.push(ob);
              fraqf_compressor_node.try_put(continue_msg{});
            }
          }
        }
      } else {
        DemuxWriterPtr w;
        auto it = fastq_writer_map.find(fname);
        if (it != fastq_writer_map.end()) {
          w = it->second;
        } else {
          int compress_level = ends_with(fname, ".gz") ? config.gzip_compress_level :
                               ends_with(fname, ".zst") ? config.zstd_compress_level :
                               1;
          w = std::make_shared<DemuxWriter>(fname, config.blocksize, compress_level);
          fastq_writer_map[fname] = w;
        }
        for (auto& r : batch) {
          size_t sz = w->add_read(std::move(r));
          if (sz >= config.blocksize) {
            w->push_current_block();
            if (!w->in_flight.exchange(true, std::memory_order_acq_rel)) {
              fastq_write_queue.push(w.get());
              fastq_writer_node.try_put(continue_msg{});
            }
          }
        }
      }
    }
    // Opportunistic local drain to prioritize writing before releasing tokens
    this->fastq_writer_node_body();
    // Chain compression -> writing in single-threaded runs
    this->fraqf_compressor_node_body();
    // After demux work is scheduled/drained, return limiter token and request next
    primary_limiter.decrement.try_put(continue_msg());
    primary_limiter.try_put(continue_msg());
    return continue_msg{};
  }

  // Write queued records for one writer
  continue_msg fastq_writer_node_body() {
    DemuxWriter* wptr = nullptr;
    while (fastq_write_queue.try_pop(wptr)) {
      if (!wptr) continue;
      std::vector<Read> block;
      while (wptr->queue.try_pop(block)) {
        wptr->writer.write_block(std::move(block));
      }
      wptr->in_flight.store(false, std::memory_order_release); // release token for this writer
    }
    return continue_msg{};
  }

  // Drain compressor queue; chain to writer queue and opportunistically drain writers
  continue_msg fraqf_compressor_node_body() {
    FraqfOutputBlockPtr ob;
    while (fraqf_compress_queue.try_pop(ob)) {
      if (!ob) continue;
      FraqfDemuxWriter* w = ob->w;
      CompressedReadBlock compressed_block = compress_reads(ob->reads, config.fraq_compress_level);
      w->add_compressed_block(ob->index, std::move(compressed_block));
      if (!w->in_flight.exchange(true, std::memory_order_acq_rel)) {
        fraqf_write_queue.push(w);
        fraqf_writer_node.try_put(continue_msg{});
      }
    }
    // Work-steal into writer to keep tail pressure
    this->fraqf_writer_node_body();
    return continue_msg{};
  }

  // Drain writer queue; each writer progresses in order using next_expected_index
  continue_msg fraqf_writer_node_body() {
    FraqfDemuxWriter* w = nullptr;
    while (fraqf_write_queue.try_pop(w)) {
      if (!w) continue;
      w->write_next_compressed_blocks();
      w->in_flight.store(false, std::memory_order_release); // release token
    }
    return continue_msg{};
  }
};

}; // namespace fraq_internal
#endif // include guard
