#ifndef _FRAQ_RUN_GRAPH_R_H
#define _FRAQ_RUN_GRAPH_R_H

#include <Rcpp.h>
#include <algorithm>
#include <iterator>
#include "fraq_run_graph.h"

namespace fraq_internal {

struct RKernelBlockInput {
  Rcpp::List reads;
  Rcpp::NumericVector index;
  size_t block_index;
};

inline RKernelBlockInput r_kernel_block_input(const BlockPtrVec &b,
                                              const FraqRunConfig &config) {
  if (b.blocks.empty() || !b.blocks[0]) {
    throw std::runtime_error("Internal error: empty R kernel input block");
  }

  const size_t block_index = b.blocks[0]->index;
  size_t nreads = b.min_reads();
  const size_t first_index = block_index * config.blocksize;
  if (config.limit > 0) {
    if (first_index >= config.limit) {
      nreads = 0;
    } else {
      nreads = std::min(nreads, config.limit - first_index);
    }
  }

  Rcpp::List reads(static_cast<R_xlen_t>(b.blocks.size()));
  Rcpp::CharacterVector read_list_names(static_cast<R_xlen_t>(b.blocks.size()));
  Rcpp::NumericVector index(static_cast<R_xlen_t>(nreads));

  for (size_t stream = 0; stream < b.blocks.size(); ++stream) {
    if (!b.blocks[stream]) {
      throw std::runtime_error("Internal error: missing R kernel input stream");
    }
    Rcpp::CharacterVector read_names(static_cast<R_xlen_t>(nreads));
    Rcpp::CharacterVector seq(static_cast<R_xlen_t>(nreads));
    Rcpp::CharacterVector qual(static_cast<R_xlen_t>(nreads));
    for (size_t i = 0; i < nreads; ++i) {
      const Read &r = b.blocks[stream]->reads[i];
      read_names[static_cast<R_xlen_t>(i)] = r.name;
      seq[static_cast<R_xlen_t>(i)] = r.seq;
      qual[static_cast<R_xlen_t>(i)] = r.qual;
      if (stream == 0) {
        index[static_cast<R_xlen_t>(i)] = static_cast<double>(first_index + i);
      }
    }
    reads[static_cast<R_xlen_t>(stream)] = Rcpp::DataFrame::create(
      Rcpp::Named("name") = read_names,
      Rcpp::Named("seq") = seq,
      Rcpp::Named("qual") = qual,
      Rcpp::Named("stringsAsFactors") = false
    );
    read_list_names[static_cast<R_xlen_t>(stream)] = "read" + std::to_string(stream + 1);
  }
  reads.names() = read_list_names;

  return RKernelBlockInput{reads, index, block_index};
}

inline Rcpp::CharacterVector r_kernel_read_column(const Rcpp::DataFrame &df,
                                                  const char *column,
                                                  const std::string &output_name) {
  if (!df.containsElementNamed(column)) {
    Rcpp::stop("fraq_run_r kernel output '%s' must contain a '%s' column",
               output_name.c_str(), column);
  }
  Rcpp::RObject value = df[column];
  if (TYPEOF(value) != STRSXP) {
    Rcpp::stop("fraq_run_r kernel output '%s' column '%s' must be character",
               output_name.c_str(), column);
  }
  return Rcpp::CharacterVector(value);
}

inline std::vector<Read> r_kernel_dataframe_to_reads(SEXP df_sexp,
                                                     const std::string &output_name) {
  if (!Rf_inherits(df_sexp, "data.frame")) {
    Rcpp::stop("fraq_run_r kernel output '%s' must be a data frame",
               output_name.c_str());
  }

  Rcpp::DataFrame df(df_sexp);
  Rcpp::CharacterVector read_names = r_kernel_read_column(df, "name", output_name);
  Rcpp::CharacterVector seq = r_kernel_read_column(df, "seq", output_name);
  Rcpp::CharacterVector qual = r_kernel_read_column(df, "qual", output_name);
  const R_xlen_t nreads = read_names.size();
  if (seq.size() != nreads || qual.size() != nreads) {
    Rcpp::stop("fraq_run_r kernel output '%s' columns must have the same length",
               output_name.c_str());
  }

  std::vector<Read> out;
  out.reserve(static_cast<size_t>(nreads));
  for (R_xlen_t i = 0; i < nreads; ++i) {
    if (Rcpp::CharacterVector::is_na(read_names[i]) ||
        Rcpp::CharacterVector::is_na(seq[i]) ||
        Rcpp::CharacterVector::is_na(qual[i])) {
      Rcpp::stop("fraq_run_r kernel output '%s' cannot contain NA reads",
                 output_name.c_str());
    }
    Read r;
    r.name = Rcpp::as<std::string>(read_names[i]);
    r.seq = Rcpp::as<std::string>(seq[i]);
    r.qual = Rcpp::as<std::string>(qual[i]);
    if (r.seq.size() != r.qual.size()) {
      Rcpp::stop("fraq_run_r kernel output '%s' has sequence and quality strings with different lengths",
                 output_name.c_str());
    }
    out.push_back(std::move(r));
  }
  return out;
}

inline ProcessedBlockPtr r_kernel_result_to_processed(SEXP result,
                                                      const size_t block_index) {
  auto p = std::make_shared<ProcessedBlock>();
  p->index = block_index;
  if (Rf_isNull(result)) {
    return p;
  }
  if (Rf_inherits(result, "data.frame") || TYPEOF(result) != VECSXP) {
    Rcpp::stop("fraq_run_r kernel must return NULL or a named list of data frames");
  }

  Rcpp::List output(result);
  SEXP output_names_sexp = Rf_getAttrib(output, R_NamesSymbol);
  if (Rf_isNull(output_names_sexp) || TYPEOF(output_names_sexp) != STRSXP) {
    Rcpp::stop("fraq_run_r kernel must return a named list of data frames");
  }
  Rcpp::CharacterVector output_names(output_names_sexp);
  if (output_names.size() != output.size()) {
    Rcpp::stop("fraq_run_r kernel must return a named list of data frames");
  }

  for (R_xlen_t i = 0; i < output.size(); ++i) {
    SEXP element = output[i];
    if (Rf_isNull(element)) {
      continue;
    }
    if (Rcpp::CharacterVector::is_na(output_names[i])) {
      Rcpp::stop("fraq_run_r kernel output names cannot be NA");
    }
    std::string output_name = Rcpp::as<std::string>(output_names[i]);
    if (output_name.empty()) {
      Rcpp::stop("fraq_run_r kernel must return a named list of data frames");
    }

    std::vector<Read> reads = r_kernel_dataframe_to_reads(element, output_name);
    if (!reads.empty()) {
      std::vector<Read> &bucket = p->bins[output_name];
      bucket.reserve(bucket.size() + reads.size());
      bucket.insert(bucket.end(),
                    std::make_move_iterator(reads.begin()),
                    std::make_move_iterator(reads.end()));
    }
  }
  return p;
}

inline ProcessedBlockPtr r_kernel_process_block(const BlockPtrVec &b,
                                                Rcpp::Function &kernel,
                                                const FraqRunConfig &config) {
  RKernelBlockInput input = r_kernel_block_input(b, config);
  Rcpp::RObject result = kernel(input.reads, input.index);
  return r_kernel_result_to_processed(result, input.block_index);
}

struct FraqRunGraphR {

  static constexpr size_t MAX_INFLIGHT_INDICES = FraqRunGraph::MAX_INFLIGHT_INDICES;

  graph flow_graph;
  const FraqRunConfig config;
  std::unordered_map<std::string, std::vector<CompressedReadBlock>> & mem_store;
  std::atomic<bool> interrupted{false};

  std::vector<MultiReader> readers;
  std::vector<size_t> current_index;
  std::vector<size_t> read_count_per_reader;
  std::atomic<size_t> block_index_limit{std::numeric_limits<size_t>::max()};

  limiter_node<continue_msg> primary_limiter;
  function_node<continue_msg, BlockPtr> primary_reader_node;
  std::vector< function_node<continue_msg, BlockPtr> > secondary_reader_nodes;

  std::unordered_map<size_t, BlockPtrVec> read_blocks;
  function_node<BlockPtr, continue_msg> joiner_node;
  sequencer_node<BlockPtrVec> kernel_queue;

  sequencer_node<ProcessedBlockPtr> reorder_node;

  std::unordered_map<std::string, DemuxWriterPtr> fastq_writer_map;
  std::unordered_map<std::string, FraqfDemuxWriterPtr> fraqf_writer_map;
  function_node<ProcessedBlockPtr, continue_msg> demux_node;
  tbb::concurrent_queue<DemuxWriter*> fastq_write_queue;

  function_node<continue_msg, continue_msg> fastq_writer_node;

  tbb::concurrent_queue<FraqfOutputBlockPtr> fraqf_compress_queue;
  tbb::concurrent_queue<FraqfDemuxWriter*> fraqf_write_queue;
  function_node<continue_msg, continue_msg> fraqf_compressor_node;
  function_node<continue_msg, continue_msg> fraqf_writer_node;

  FraqRunGraphR(const std::vector<std::string> &input_files,
                const FraqRunConfig config,
                std::unordered_map<std::string, std::vector<CompressedReadBlock>> & mem_store) :
    config(config),
    mem_store(mem_store),
    readers(),
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
    kernel_queue(flow_graph,
      [](const BlockPtrVec &b) -> size_t {
        return b.blocks[0]->index;
      }
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
    readers = readers_init(input_files);
    const size_t limit_block = (config.limit > 0)
      ? ((config.limit - 1) / config.blocksize)
      : std::numeric_limits<size_t>::max();
    block_index_limit.store(limit_block, std::memory_order_release);

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

    make_edge(primary_limiter, primary_reader_node);
    make_edge(primary_reader_node, joiner_node);
    for (auto &fn : secondary_reader_nodes) {
      make_edge(fn, joiner_node);
    }
    make_edge(reorder_node, demux_node);

    for(size_t k = 0; k < MAX_INFLIGHT_INDICES; ++k) {
      primary_limiter.try_put(continue_msg());
    }
  }

  bool try_pop_kernel_block(BlockPtrVec &b) {
    return kernel_queue.try_get(b);
  }

  void submit_processed_block(const ProcessedBlockPtr &p) {
    if (!p) {
      throw std::runtime_error("Internal error: null processed block");
    }
    if (!reorder_node.try_put(p)) {
      throw std::runtime_error("Internal error: failed to enqueue processed block");
    }
  }

  void wait_for_idle() {
    flow_graph.wait_for_all();
  }

  void request_interrupt() {
    interrupted.store(true, std::memory_order_release);
  }

  void wait_and_flush() {
    flow_graph.wait_for_all();
    if (was_interrupted()) {
      cleanup_outputs();
      return;
    }

    for (auto& kv : fraqf_writer_map) {
      FraqfDemuxWriter* w = kv.second.get();
      if (w->current_block_size() > 0) {
        FraqfOutputBlockPtr ob = w->release_current_block();
        fraqf_compress_queue.push(ob);
        fraqf_compressor_node.try_put(continue_msg{});
      }
    }
    this->fraqf_compressor_node_body();

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
    flow_graph.wait_for_all();

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

  bool was_interrupted() const {
    return interrupted.load(std::memory_order_acquire);
  }

  void cleanup_outputs() {
    for (auto& kv : fastq_writer_map) {
      const std::string& path = kv.first;
      if (ends_with(path, ".fifo")) {
        try {
          kv.second->writer.close();
        } catch (...) {
        }
        continue;
      }
      try {
        kv.second->writer.close();
      } catch (...) {
      }
      std::remove(path.c_str());
    }
    for (auto& kv : fraqf_writer_map) {
      const std::string& path = kv.first;
      if (ends_with(path, ".mem")) {
        mem_store.erase(path);
        continue;
      }
      try {
        kv.second->writer.close();
      } catch (...) {
      }
      std::remove(path.c_str());
    }
  }

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

  bool interrupt_requested() {
    if (interrupted.load(std::memory_order_acquire)) {
      return true;
    }
    if (config.interrupt_fn && config.interrupt_fn(config.interrupt_ctx)) {
      interrupted.store(true, std::memory_order_release);
      return true;
    }
    return false;
  }

  BlockPtr primary_reader_step() {
      if (interrupt_requested()) {
        auto b = std::make_shared<Block>();
        b->index = current_index[0];
        b->reader_index = 0;
        return b;
      }
      for (auto &fn : secondary_reader_nodes) {
        fn.try_put(continue_msg{});
      }
      return reader_node_body(0);
  }

  BlockPtr reader_node_body(const size_t ri) {
      size_t& idx = current_index[ri];
      const size_t limit = block_index_limit.load(std::memory_order_acquire);
      auto b = std::make_shared<Block>();
      b->index = idx;
      b->reader_index = ri;
      if (idx > limit) {
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
    if (!b || b->reads.empty()) {
      return continue_msg();
    }
    size_t idx = b->index;
    auto it = read_blocks.find(idx);
    if (it == read_blocks.end()) {
      BlockPtrVec bpv;
      bpv.resize(readers.size());
      auto ins = read_blocks.emplace(idx, std::move(bpv));
      it = ins.first;
    }
    it->second.insert(b, b->reader_index);

    if(it->second.full()) {
      BlockPtrVec to_process = std::move(it->second);
      read_blocks.erase(it);
      if (to_process.blocks.empty() || !to_process.blocks[0]) {
        throw std::runtime_error("Internal error: invalid R kernel block sequence");
      }
      if (!kernel_queue.try_put(to_process)) {
        throw std::runtime_error("Internal error: failed to enqueue R kernel block");
      }
    }
    return continue_msg();
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
    this->fastq_writer_node_body();
    this->fraqf_compressor_node_body();
    tbb_compat::decrementer(primary_limiter).try_put(continue_msg());
    primary_limiter.try_put(continue_msg());
    return continue_msg{};
  }

  continue_msg fastq_writer_node_body() {
    DemuxWriter* wptr = nullptr;
    while (fastq_write_queue.try_pop(wptr)) {
      if (!wptr) continue;
      std::vector<Read> block;
      while (wptr->queue.try_pop(block)) {
        wptr->writer.write_block(std::move(block));
      }
      wptr->in_flight.store(false, std::memory_order_release);
    }
    return continue_msg{};
  }

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
    this->fraqf_writer_node_body();
    return continue_msg{};
  }

  continue_msg fraqf_writer_node_body() {
    FraqfDemuxWriter* w = nullptr;
    while (fraqf_write_queue.try_pop(w)) {
      if (!w) continue;
      w->write_next_compressed_blocks();
      w->in_flight.store(false, std::memory_order_release);
    }
    return continue_msg{};
  }
};

}; // namespace fraq_internal
#endif // include guard
