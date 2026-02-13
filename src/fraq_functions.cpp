#include <Rcpp.h>

#include "fraq_defines.h"
#include "fraq_run_graph.h"
#include "fraq_concat.h"
#include "distance_functions.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <atomic>
#include <chrono>
#include <exception>
#include <thread>
#include <unordered_set>
#include <tbb/global_control.h>
#include <tbb/enumerable_thread_specific.h>

using namespace Rcpp;
using fraq::output_t; // std::vector<outpair_t>, outpair_t = std::pair<std::string /*output file*/, Read>
using fraq::input_t; // std::vector<fraq::Read>
using fraq::zip; // zip(vector<A>, vector<B>) -> vector<pair<A,B>>
using fraq::ends_with;

// Forward declarations for exported functions
void fraq_run(const std::vector<std::string> &input_files, fraq::process_task_t task, const int nthreads, fraq::FraqRunConfig config);
void fraq_concat(const std::vector<std::string> &input_files, const std::string &output_file, const int nthreads, const fraq::FraqRunConfig config);
int fraq_options(const std::string &option, const int value, const bool set);
void fraq_export_functions(DllInfo *dll);
int rcpp_fraq_options(const std::string &option, const int value, const bool set);
void rcpp_fraq_downsample(std::vector<std::string> input, const std::vector<std::string> &output, const double amount, const int nthreads);
void rcpp_fraq_slice(std::vector<std::string> input, const std::vector<std::string> &output, const double limit, const std::vector<double> &select, const int nthreads);
void rcpp_fraq_convert(std::vector<std::string> input, const std::vector<std::string> &output, const int nthreads);
void rcpp_fraq_chunk(std::vector<std::string> input, const std::vector<std::string> &output_prefix, const std::string &output_suffix, const double chunk_size, const int nthreads);
void rcpp_fraq_concat(const std::vector<std::string> &input, const std::string &output, const int nthreads);
DataFrame rcpp_fraq_count_barcodes(const std::vector<std::string> &input, const std::vector<std::string> &barcodes, const int max_distance, const bool allow_revcomp, const int nthreads);
void rcpp_fraq_demux(const std::vector<std::string> &input, const std::vector<std::string> &output_format, const std::vector<std::string> &barcodes, const int max_distance, const int nthreads);
Rcpp::List rcpp_fraq_merge_pairs(const std::vector<std::string> &input, const std::string &output_merged, const std::vector<std::string> &output_unmerged, const int min_overlap, const double max_mismatch_rate, const std::string &consensus_mode, const bool trim_overhang, const bool revcomp_R2, const int nthreads);
void rcpp_fraq_quality_filter(const std::vector<std::string> &input, const std::vector<std::string> &output, const double min_mean_quality, const int max_low_q_bases, const int low_q_threshold, const int nthreads);
DataFrame rcpp_fraq_trim_adapters(const std::vector<std::string> &input, const std::vector<std::string> &output, const std::vector<std::string> &adapters, const int max_distance, const bool filter_untrimmed, const int nthreads);
Rcpp::List rcpp_fraq_summary(const std::vector<std::string> &input, const bool phred33, const int min_overlap, const double max_mismatch_rate, const size_t limit, const int nthreads);
void rcpp_fraq_wait(std::vector<std::string> input, const std::vector<std::string> &output, const int sleep_ms, const int nthreads);
DataFrame rcpp_fraq_align(CharacterVector query, CharacterVector target, int max_distance, std::string ambiguity_base, std::string boundary, std::string distance_metric);
DataFrame rcpp_fraq_mem_list();
bool rcpp_fraq_mem_remove(const std::string &mem_key);

// global variables, set using fraq_options()
// max allowed value of 2^16-1 = 65535
static size_t BLOCKSIZE = 65535;
static int FRAQ_COMPRESS_LEVEL = 3;
static int ZSTD_COMPRESS_LEVEL = 3;
static int GZIP_COMPRESS_LEVEL = 6;
// Note: kernel parallelism is now chosen per-operation; no global serial flag.

// global in memory fastq store
static std::unordered_map<std::string, std::vector<fraq_internal::CompressedReadBlock>> mem_store;

struct InterruptState {
  std::atomic<bool> requested{false};
  std::exception_ptr exception;
  std::thread::id main_thread_id;
};

static bool interrupt_r_check(void* ctx) {
  auto* state = static_cast<InterruptState*>(ctx);
  if (!state) return false;
  if (state->main_thread_id != std::this_thread::get_id()) {
    return false;
  }
  try {
    Rcpp::checkUserInterrupt();
  } catch (...) {
    state->requested.store(true, std::memory_order_release);
    state->exception = std::current_exception();
    return true;
  }
  return false;
}

// Helper used by multiple kernels to compute the best overlap length
static size_t fraq_calculate_overlap(const std::string &seq1,
                                     const std::string &seq2,
                                     int min_overlap,
                                     double max_mismatch_rate,
                                     int &mismatches,
                                     bool treat_n_as_mismatch = true) {
  const int L1 = static_cast<int>(seq1.size());
  const int L2 = static_cast<int>(seq2.size());
  const int maxL = std::min(L1, L2);
  bool too_many = false;

  for (int L = maxL; L >= min_overlap; --L) {
    int start1 = L1 - L;
    int mismatch_here = 0;
    int max_mism_allowed = static_cast<int>(std::ceil(max_mismatch_rate * static_cast<double>(L)));
    for (int k = 0; k < L; ++k) {
      char a = seq1[static_cast<size_t>(start1 + k)];
      char b = seq2[static_cast<size_t>(k)];
      bool is_N = treat_n_as_mismatch && (a == 'N' || a == 'n' || b == 'N' || b == 'n');
      if (is_N || a != b) {
        ++mismatch_here;
        if (mismatch_here > max_mism_allowed) {
          mismatch_here = -1;
          too_many = true;
          break;
        }
      }
    }
    if (mismatch_here >= 0) {
      mismatches = mismatch_here;
      return static_cast<size_t>(L);
    }
  }
  mismatches = too_many ? -1 : 0;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// functions exported to C++

void fraq_run(const std::vector<std::string> &input_files, fraq::process_task_t task, int nthreads = 1, fraq::FraqRunConfig config = fraq::FraqRunConfig{}) {
  if (nthreads <= 1) {
    nthreads = 1;
  }
  InterruptState interrupt_state;
  interrupt_state.main_thread_id = std::this_thread::get_id();
  config.interrupt_fn = interrupt_r_check;
  config.interrupt_ctx = &interrupt_state;
  tbb::global_control gc(tbb::global_control::parameter::max_allowed_parallelism, nthreads);
  fraq_internal::FraqRunGraph fg(input_files, task, config, mem_store);
  fg.wait_and_flush();
  if (interrupt_state.exception) {
    Rcpp::warning("fraq interrupted");
    std::rethrow_exception(interrupt_state.exception);
  } else {
    auto excess = fg.check_reader_balance();
    if (!excess.empty()) {
      Rcpp::Rcerr << "fraq_run detected excess reads from inputs: ";
      bool first = true;
      for (size_t idx : excess) {
        if (!first) Rcpp::Rcerr << ", ";
        Rcpp::Rcerr << idx;
        first = false;
      }
      Rcpp::Rcerr << std::endl;
    }
  }
}

inline void fraq_concat(const std::vector<std::string> &input_files,
                        const std::string &output_file,
                        const int nthreads,
                        const fraq::FraqRunConfig config = fraq::FraqRunConfig{}) {
  if (input_files.empty()) {
    throw std::runtime_error("fraq_concat requires at least one input file");
  }
  if (output_file.empty()) {
    throw std::runtime_error("fraq_concat output path must be non-empty");
  }
  tbb::global_control gc(tbb::global_control::parameter::max_allowed_parallelism, nthreads);
  const bool is_fraq_output = ends_with(output_file, ".fraq") || ends_with(output_file, ".mem");
  if (is_fraq_output) {
    fraq_internal::FraqConcatGraph fg(input_files, config, mem_store, output_file);
    fg.start();
    fg.wait();
    fg.flush();
  } else {
    fraq_internal::FastqConcatGraph fg(input_files, config, mem_store, output_file);
    fg.start();
    fg.wait();
    fg.flush();
  }
}


int fraq_options(const std::string & option, const int value, const bool set) {
  if(option == "blocksize") {
    int output = static_cast<int>(BLOCKSIZE);
    if(set) {
      if(value < 1 || value > 65535) throw(std::runtime_error("blocksize must be between 1 and 65535"));
      BLOCKSIZE = static_cast<size_t>(value);
    }
    return output;
  } else if(option == "fraq_compress_level") {
    int output = FRAQ_COMPRESS_LEVEL;
    if(set) {
      if(value < 1 || value > 22) throw(std::runtime_error("fraq_compress_level must be between 1 and 22"));
      FRAQ_COMPRESS_LEVEL = value;
    }
    return output;
  } else if(option == "zstd_compress_level") {
    int output = ZSTD_COMPRESS_LEVEL;
    if(set) {
      if(value < 1 || value > 22) throw(std::runtime_error("zstd_compress_level must be between 1 and 22"));
      ZSTD_COMPRESS_LEVEL = value;
    }
    return output;
  } else if(option == "gzip_compress_level") {
    int output = GZIP_COMPRESS_LEVEL;
    if(set) {
      if(value < 1 || value > 9) throw(std::runtime_error("gzip_compress_level must be between 1 and 9"));
      GZIP_COMPRESS_LEVEL = value;
    }
    return output;
  } else {
    throw(std::runtime_error("unknown option, valid options are blocksize, fraq_compress_level, zstd_compress_level, gzip_compress_level"));
  }
}

// [[Rcpp::init]]
void fraq_export_functions(DllInfo* dll) {
    R_RegisterCCallable("fraq", "fraq_run",
    (DL_FUNC)static_cast<void(*)(const std::vector<std::string>&,
                                  fraq::process_task_t,
                                  const int,
                                  fraq::FraqRunConfig)> (&fraq_run) );
    R_RegisterCCallable("fraq", "fraq_concat",
    (DL_FUNC)static_cast<void(*)(const std::vector<std::string>&,
                                  const std::string&,
                                  const int,
                                  const fraq::FraqRunConfig)> (&fraq_concat) );
    R_RegisterCCallable("fraq", "fraq_options", (DL_FUNC)&fraq_options);
    R_RegisterCCallable("fraq", "fraq_hm_starts", (DL_FUNC)&fraq_hm_starts);
    R_RegisterCCallable("fraq", "fraq_hm_contains", (DL_FUNC)&fraq_hm_contains);
    R_RegisterCCallable("fraq", "fraq_hm_global", (DL_FUNC)&fraq_hm_global);
    R_RegisterCCallable("fraq", "fraq_lv_starts", (DL_FUNC)&fraq_lv_starts);
    R_RegisterCCallable("fraq", "fraq_lv_contains", (DL_FUNC)&fraq_lv_contains);
    R_RegisterCCallable("fraq", "fraq_lv_global", (DL_FUNC)&fraq_lv_global);
}

///////////////////////////////////////////////////////////////////////////////
// functions exported to R package

// [[Rcpp::export(rng=false)]]
int rcpp_fraq_options(const std::string & option, const int value, const bool set) {
  return fraq_options(option, value, set);
}

// [[Rcpp::export(rng=false)]]
void rcpp_fraq_downsample(std::vector<std::string> input, const std::vector<std::string> &output, const double amount, const int nthreads = 1) {
  if(input.size() == 0 || input.size() != output.size()) throw(std::runtime_error("input and output must be character vectors of the same length > 0"));
  if(amount > 1 || amount < 0) throw(std::runtime_error("amount must between 0 and 1"));
  auto downsample_task = [&](input_t reads, size_t index) -> output_t {
    const double scaled_index = amount * static_cast<double>(index);
    const double scaled_next = amount * static_cast<double>(index + 1);
    if (std::floor(scaled_next) > std::floor(scaled_index)) {
      return zip(output, std::move(reads));
    }
    return {};
  };
  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_run(input, downsample_task, nthreads, config);
}

// [[Rcpp::export(rng=false)]]
void rcpp_fraq_slice(std::vector<std::string> input,
                     const std::vector<std::string> &output,
                     const size_t limit,
                     const std::vector<size_t> &select,
                     const int nthreads = 1) {
  if (input.empty() || input.size() != output.size()) {
    throw std::runtime_error("input and output must be character vectors of the same length > 0");
  }

  bool has_select = limit == 0; // assume we want to use select if limit not supplied (0 is a special value that indicates no limit)
  if (has_select) {
    std::unordered_set<size_t> select_ids;
    if (select.size() == 0) { // early exit for special case where numeric(0) is supplied to select
      return;
    }
    select_ids.reserve(select.size());
    size_t max_id = 0;
    for (const size_t id : select) {
      select_ids.insert(id);
      if (id > max_id) {
        max_id = id;
      }
    }
    auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false, max_id + 1);
    auto select_task = [&](input_t reads, size_t index) -> output_t {
      const bool keep = select_ids.find(index) != select_ids.end();
      if (!keep) {
        return {};
      }
      return zip(output, std::move(reads));
    };
    fraq_run(input, select_task, nthreads, config);
  } else {
    auto identity_task = [&](input_t reads, size_t index) -> output_t {
      return zip(output, std::move(reads));
    };
    auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false, limit);
    fraq_run(input, identity_task, nthreads, config);
  }
}

// [[Rcpp::export(rng=false)]]
void rcpp_fraq_convert(std::vector<std::string> input, const std::vector<std::string> &output, const int nthreads = 1) {
  if(input.size() == 0 || input.size() != output.size()) throw(std::runtime_error("input and output must be character vectors of the same length > 0"));
  auto convert_task = [&](input_t reads, size_t index) -> output_t {
    return zip(output, std::move(reads));
  };
  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_run(input, convert_task, nthreads, config);
}

// [[Rcpp::export(rng=false)]]
void rcpp_fraq_chunk(std::vector<std::string> input, const std::vector<std::string> &output_prefix, const std::string &output_suffix, const double chunk_size, const int nthreads = 1) {
  size_t chunk_size_value = static_cast<size_t>(chunk_size);
  if(chunk_size <= 0 || chunk_size_value == 0) throw(std::runtime_error("chunk_size must be >= 1"));
  if(input.size() == 0 || output_prefix.size() == 0) throw(std::runtime_error("input and output must be character vectors of the same length > 0"));
  if(input.size() != output_prefix.size()) throw(std::runtime_error("input and output must be character vectors of the same length > 0"));

  const std::string resolved_suffix = [&]() -> std::string {
    if (output_suffix == "zst") return ".fastq.zst";
    if (output_suffix == "gz") return ".fastq.gz";
    if (output_suffix == "fastq") return ".fastq";
    if (output_suffix == "fraq") return ".fraq";
    if (output_suffix == "mem") return ".mem";
    throw std::runtime_error("output_suffix must be one of (zst, gz, fastq, fraq or mem)");
  }();

  auto chunk_task = [&](input_t reads, size_t index) -> output_t {
    size_t chunk_number = index / chunk_size_value;
    output_t outputs(reads.size());
    outputs.reserve(output_prefix.size());
    for(size_t i = 0; i < output_prefix.size(); ++i) {
      std::string target = output_prefix[i];
      target += "_chunk";
      target += std::to_string(chunk_number);
      target += resolved_suffix;
      outputs[i] = std::make_pair( std::move(target), std::move(reads[i]) );
    }
    return outputs;
  };
  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_run(input, chunk_task, nthreads, config);
}


// [[Rcpp::export(rng=false)]]
void rcpp_fraq_concat(const std::vector<std::string> &input,
                      const std::string &output,
                      const int nthreads = 1) {
  if (output.empty()) throw std::runtime_error("output must be a non-empty string");
  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_concat(input, output, nthreads, config);
}




// [[Rcpp::export(rng=false)]]
DataFrame rcpp_fraq_count_barcodes(const std::vector<std::string> &input,
                                   const std::vector<std::string> &barcodes,
                                   const int max_distance = 1,
                                   const bool allow_revcomp = false,
                                   const int nthreads = 1) {
  if(input.size() == 0) throw(std::runtime_error("input must be a character vector of length > 0"));
  if(barcodes.size() == 0) throw(std::runtime_error("barcodes must be a character vector of length > 0"));
  if(max_distance < 0) throw(std::runtime_error("max_distance must be non-negative"));

  // per thread barcode counts
  tbb::enumerable_thread_specific<std::unordered_map<std::string, uint64_t>> count_maps;

  std::vector<std::string> barcodes_revcomp(barcodes.size());
  std::transform(barcodes.begin(), barcodes.end(), barcodes_revcomp.begin(), fraq::reverse_complement);

  // helper function
  auto do_match = [&](const std::string & seq, const std::string & bc) {
    return bc.size() <= seq.size() && fraq_hm_contains(bc, seq, max_distance).distance <= max_distance;
  };
  auto count_barcodes_task = [&](input_t reads, size_t index) -> output_t {
    std::unordered_map<std::string, uint64_t> &local_map = count_maps.local();
    std::string bc_match;
    for(size_t i= 0; i < barcodes.size(); ++i) {
      const std::string &barcode = barcodes[i];
      const std::string &barcode_revcomp = barcodes_revcomp[i];
      bool matched = false;
      for(auto & read : reads) {
        if (do_match(read.seq, barcode) ||
            (allow_revcomp && do_match(read.seq, barcode_revcomp))) {
          matched = true;
          break;
        }
      }
      if (!matched) continue;
      if (bc_match.empty()) {
        bc_match = barcode; // first match found
      } else if (bc_match != barcode) {
        local_map["MULTI_MATCH"]++;
        return {};
      }
    }
    if(!bc_match.empty()) {
      local_map[bc_match]++;
    } else {
      local_map["NO_MATCH"]++;
    } 
    return {};
  };
  {
    auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
    fraq_run(input, count_barcodes_task, nthreads, config);
  }

  // Combine results from all threads
  std::unordered_map<std::string, uint64_t> total_counts;
  for (const auto &local_map : count_maps) {
    for (const auto &pair : local_map) {
      total_counts[pair.first] += pair.second;
    }
  }

  // convert to DataFrame
  CharacterVector barcodes_out(total_counts.size());
  NumericVector counts_out(total_counts.size());
  int i = 0;
  for (const auto &pair : total_counts) {
    barcodes_out[i] = pair.first;
    counts_out[i] = static_cast<double>(pair.second);
    i++;
  }
  return DataFrame::create(
      _["barcode"] = barcodes_out,
      _["count"]   = counts_out
  );
}

// [[Rcpp::export(rng=false)]]
void rcpp_fraq_demux(const std::vector<std::string> &input,
                     const std::vector<std::string> &output_format,
                     const std::vector<std::string> &barcodes,
                     const int max_distance = 1,
                     const int nthreads = 1) {
  if(input.size() == 0) throw(std::runtime_error("input must be a character vector of length > 0"));
  if(input.size() != output_format.size()) throw(std::runtime_error("input and output_format must be the same length"));
  if(barcodes.size() == 0) throw(std::runtime_error("barcodes must be a character vector of length > 0"));
  if(max_distance < 0) throw(std::runtime_error("max_distance must be non-negative"));

  const std::string placeholder = "{barcode}";

  auto starts_match = [&](const std::string &barcode, const std::string &sequence) -> bool {
    return barcode.size() <= sequence.size() && fraq_hm_starts(barcode, sequence, max_distance).distance <= max_distance;
  };

  std::vector<std::string> output_prefix(input.size());
  std::vector<std::string> output_suffix(input.size());
  for(size_t i=0; i<input.size(); ++i) {
    const std::string &fmt = output_format[i];
    std::size_t pos = fmt.find(placeholder);
    if (pos == std::string::npos) {
      throw std::runtime_error("output format must be in the form <file_prefix>{barcode}<file_suffix>, e.g. Sample1_R1_{barcode}.fastq.zst");
    }
    output_prefix[i] = fmt.substr(0, pos);
    output_suffix[i] = fmt.substr(pos + placeholder.size());
  }

  auto zip_demux_output = [&](input_t &&reads, const std::string &barcode) {
    output_t out;
    out.reserve(output_prefix.size());
    if (reads.size() != output_prefix.size()) {
      throw std::runtime_error("fraq_demux expected reads.size() == output_format.size()");
    }
    for(size_t j = 0; j < output_prefix.size(); ++j) {
      out.emplace_back(output_prefix[j] + barcode + output_suffix[j], std::move(reads[j]));
    }
    return out;
  };
  auto demux_task = [&](input_t reads, size_t /*index*/) -> output_t {
    if (reads.empty()) return {};
    std::string bc_match;
    for(size_t i = 0; i < barcodes.size(); ++i) {
      const std::string &barcode = barcodes[i];
      // only look at the start of first read
      bool is_match = starts_match(barcode, reads[0].seq);
      if(is_match) {
        if(bc_match.empty()) {
          bc_match = barcode;
        } else {
          return zip_demux_output(std::move(reads), "MULTI_MATCH");
        }
      }
    }
    if(!bc_match.empty()) {
      return zip_demux_output(std::move(reads), bc_match);
    } else {
      return zip_demux_output(std::move(reads), "NO_MATCH");
    }
  };
  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_run(input, demux_task, nthreads, config);
}

// [[Rcpp::export(rng=false)]]
void rcpp_fraq_quality_filter(const std::vector<std::string> &input,
                              const std::vector<std::string> &output,
                              const double min_mean_quality,
                              const int max_low_q_bases,
                              const int low_q_threshold,
                              const int nthreads = 1) {
  if (input.empty()) throw std::runtime_error("input must be a character vector of length > 0");
  if (input.size() != output.size()) throw std::runtime_error("input and output must be character vectors of the same length");
  if (min_mean_quality < 0.0) throw std::runtime_error("min_mean_quality must be non-negative");
  if (max_low_q_bases < 0) throw std::runtime_error("max_low_q_bases must be non-negative");
  if (low_q_threshold < 0) throw std::runtime_error("low_q_threshold must be non-negative");

  auto passes_filters = [&](const fraq::Read &read) -> bool {
    const std::string &qual = read.qual;
    const size_t qlen = qual.size();
    if (qlen == 0) return false;

    double total_q = 0.0;
    int low_q_count = 0;
    for (char c : qual) {
      int q = static_cast<int>(static_cast<unsigned char>(c)) - 33; // PHRED+33
      if (q < 0) q = 0;
      total_q += static_cast<double>(q);
      if (q < low_q_threshold) {
        ++low_q_count;
        if (low_q_count > max_low_q_bases) {
          return false;
        }
      }
    }
    double mean_q = total_q / static_cast<double>(qlen);
    if (mean_q < min_mean_quality) return false;
    return true;
  };

  auto filter_task = [&](input_t reads, size_t /*index*/) -> output_t {
    if (reads.size() != output.size()) {
      throw std::runtime_error("fraq_quality_filter expected reads.size() == output.size()");
    }
    for (const auto &read : reads) {
      if (!passes_filters(read)) {
        return {};
      }
    }
    return zip(output, std::move(reads));
  };

  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_run(input, filter_task, nthreads, config);
}

// [[Rcpp::export(rng=false)]]
Rcpp::List rcpp_fraq_merge_pairs(const std::vector<std::string> &input,
                                 const std::string &output_merged,
                                 const std::vector<std::string> &output_unmerged,
                                 const int min_overlap = 12,
                                 const double max_mismatch_rate = 0.10,
                                 const std::string &consensus_mode = "max",
                                 const bool trim_overhang = true,
                                 const bool revcomp_R2 = true,
                                 const int nthreads = 1) {
  if (input.size() != 2) throw std::runtime_error("input must be a character vector of length 2 (R1, R2)");
  if (output_merged.empty()) throw std::runtime_error("output_merged must be a non-empty string");
  if (!(output_unmerged.empty() || output_unmerged.size() == 2)) {
    throw std::runtime_error("output_unmerged must be length 0 or length 2");
  }
  if (min_overlap < 1) throw std::runtime_error("min_overlap must be >= 1");
  if (max_mismatch_rate < 0.0 || max_mismatch_rate > 1.0) {
    throw std::runtime_error("max_mismatch_rate must be between 0 and 1");
  }
  std::string mode_lc = consensus_mode;
  std::transform(mode_lc.begin(), mode_lc.end(), mode_lc.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  enum class ConsensusMode { MAX, MEAN, R1, R2 };
  ConsensusMode mode;
  if (mode_lc == "max") {
    mode = ConsensusMode::MAX;
  } else if (mode_lc == "mean") {
    mode = ConsensusMode::MEAN;
  } else if (mode_lc == "r1") {
    mode = ConsensusMode::R1;
  } else if (mode_lc == "r2") {
    mode = ConsensusMode::R2;
  } else {
    throw std::runtime_error("consensus_mode must be one of ('max','mean','r1','r2')");
  }

  struct MergeStats {
    uint64_t merged = 0;
    uint64_t unmerged = 0;
    double sum_insert = 0.0;
    double sumsq_insert = 0.0;
    double sum_overlap = 0.0;
    double sum_mismatch_rate = 0.0;
  };

  tbb::enumerable_thread_specific<MergeStats> stats_tls;

  auto reverse_quality = [](const std::string &qual) {
    std::string out = qual;
    std::reverse(out.begin(), out.end());
    return out;
  };

  auto merge_task = [&](input_t reads, size_t /*index*/) -> output_t {
    if (reads.size() != 2) {
      throw std::runtime_error("fraq_merge_pairs expects exactly 2 reads per record");
    }
    MergeStats &local = stats_tls.local();
    const fraq::Read &r1 = reads[0];
    const fraq::Read &r2 = reads[1];

    std::string seq2 = revcomp_R2 ? fraq::reverse_complement(r2.seq) : r2.seq;
    std::string qual2 = revcomp_R2 ? reverse_quality(r2.qual) : r2.qual;

    int mismatches = 0;
    size_t overlap = fraq_calculate_overlap(r1.seq, seq2, min_overlap, max_mismatch_rate, mismatches, true);
    if (overlap == 0) {
      local.unmerged++;
      if (output_unmerged.empty()) {
        return {};
      }
      return zip(output_unmerged, std::move(reads));
    }

    const size_t start1 = r1.seq.size() - static_cast<size_t>(overlap);
    const size_t overlap_end1 = start1 + static_cast<size_t>(overlap);

    std::string merged_seq;
    std::string merged_qual;
    if (trim_overhang) {
      merged_seq.reserve(start1 + static_cast<size_t>(overlap) + (seq2.size() > static_cast<size_t>(overlap) ? seq2.size() - static_cast<size_t>(overlap) : 0));
      merged_qual.reserve(merged_seq.capacity());
      merged_seq.append(r1.seq.begin(), r1.seq.begin() + static_cast<std::ptrdiff_t>(start1));
      merged_qual.append(r1.qual.begin(), r1.qual.begin() + static_cast<std::ptrdiff_t>(start1));
    } else {
      merged_seq.reserve(static_cast<size_t>(overlap));
      merged_qual.reserve(static_cast<size_t>(overlap));
    }

    for (size_t k = 0; k < overlap; ++k) {
      size_t idx1 = start1 + k;
      char b1 = r1.seq[idx1];
      char b2 = seq2[k];
      char q1 = r1.qual[idx1];
      char q2 = qual2[k];

      char chosen_base;
      char chosen_q;

      switch (mode) {
        case ConsensusMode::R1:
          chosen_base = b1;
          chosen_q = q1;
          break;
        case ConsensusMode::R2:
          chosen_base = b2;
          chosen_q = q2;
          break;
        case ConsensusMode::MEAN:
          if (b1 == b2) {
            chosen_base = b1;
          } else if (q1 >= q2) {
            chosen_base = b1;
          } else {
            chosen_base = b2;
          }
          chosen_q = static_cast<char>(std::round((static_cast<double>(static_cast<unsigned char>(q1)) + static_cast<double>(static_cast<unsigned char>(q2))) / 2.0));
          break;
        case ConsensusMode::MAX:
        default:
          if (b1 == b2) {
            chosen_base = b1;
            chosen_q = (q1 >= q2) ? q1 : q2;
          } else if (q1 >= q2) {
            chosen_base = b1;
            chosen_q = q1;
          } else {
            chosen_base = b2;
            chosen_q = q2;
          }
          break;
      }

      merged_seq.push_back(chosen_base);
      merged_qual.push_back(chosen_q);
    }

    if (trim_overhang) {
      if (overlap_end1 < r1.seq.size()) {
        merged_seq.append(r1.seq.begin() + static_cast<std::ptrdiff_t>(overlap_end1), r1.seq.end());
        merged_qual.append(r1.qual.begin() + static_cast<std::ptrdiff_t>(overlap_end1), r1.qual.end());
      }
      if (static_cast<size_t>(overlap) < seq2.size()) {
        merged_seq.append(seq2.begin() + static_cast<std::ptrdiff_t>(overlap), seq2.end());
        merged_qual.append(qual2.begin() + static_cast<std::ptrdiff_t>(overlap), qual2.end());
      }
    }

    fraq::Read merged_read = r1;
    merged_read.seq = std::move(merged_seq);
    merged_read.qual = std::move(merged_qual);

    double insert_size = static_cast<double>(r1.seq.size() + r2.seq.size() - static_cast<size_t>(overlap));
    local.merged++;
    local.sum_insert += insert_size;
    local.sumsq_insert += insert_size * insert_size;
    local.sum_overlap += static_cast<double>(overlap);
    local.sum_mismatch_rate += overlap ? static_cast<double>(mismatches) / static_cast<double>(overlap) : 0.0;

    output_t out;
    out.reserve(1);
    out.emplace_back(output_merged, std::move(merged_read));
    return out;
  };

  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_run(input, merge_task, nthreads, config);

  MergeStats total_stats;
  for (const auto &s : stats_tls) {
    total_stats.merged += s.merged;
    total_stats.unmerged += s.unmerged;
    total_stats.sum_insert += s.sum_insert;
    total_stats.sumsq_insert += s.sumsq_insert;
    total_stats.sum_overlap += s.sum_overlap;
    total_stats.sum_mismatch_rate += s.sum_mismatch_rate;
  }

  double mean_insert = total_stats.merged ? (total_stats.sum_insert / static_cast<double>(total_stats.merged)) : NA_REAL;
  double sd_insert = NA_REAL;
  if (total_stats.merged > 1) {
    double mean_sq = total_stats.sumsq_insert / static_cast<double>(total_stats.merged);
    double var = mean_sq - mean_insert * mean_insert;
    if (var < 0.0) var = 0.0;
    sd_insert = std::sqrt(var);
  }

  double mean_overlap = total_stats.merged ? (total_stats.sum_overlap / static_cast<double>(total_stats.merged)) : NA_REAL;
  double mean_mismatch_rate = total_stats.merged ? (total_stats.sum_mismatch_rate / static_cast<double>(total_stats.merged)) : NA_REAL;

  return Rcpp::List::create(
      _["merged_reads"] = static_cast<double>(total_stats.merged),
      _["unmerged_reads"] = static_cast<double>(total_stats.unmerged),
      _["mean_insert_size"] = mean_insert,
      _["sd_insert_size"] = sd_insert,
      _["mean_overlap"] = mean_overlap,
      _["mean_mismatch_rate"] = mean_mismatch_rate
  );
}

// [[Rcpp::export(rng=false)]]
DataFrame rcpp_fraq_trim_adapters(const std::vector<std::string> &input,
                                  const std::vector<std::string> &output,
                                  const std::vector<std::string> &adapters,
                                  const int max_distance = 1,
                                  const bool filter_untrimmed = true,
                                  const int nthreads = 1) {
  if(input.size() == 0 || input.size() != output.size()) throw(std::runtime_error("input and output must be character vectors of the same length > 0"));
  if(adapters.size() == 0) throw(std::runtime_error("adapters must be a character vector of length > 0"));
  if(max_distance < 0) throw(std::runtime_error("max_distance must be non-negative"));

  tbb::enumerable_thread_specific<std::unordered_map<std::string, uint64_t>> count_maps;
  auto trim_task = [&](input_t reads, size_t index) -> output_t {
    std::unordered_map<std::string, uint64_t> &local_map = count_maps.local();
    for(auto & adapter : adapters) {
      auto res = fraq_hm_starts(adapter, reads[0].seq, max_distance);
      if(res.distance <= max_distance) {
        const size_t trim_start = static_cast<size_t>(res.end + 1);
        const size_t trim_end = reads[0].seq.size();
        fraq::trim_read(reads[0], trim_start, trim_end);
        local_map[adapter]++;
        return zip(output, std::move(reads));
      }
    }
    local_map["NO_ADAPTER"]++;
    if(filter_untrimmed) {
      return {};
    } else {
      return zip(output, std::move(reads));
    }
  };
  {
    auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
    fraq_run(input, trim_task, nthreads, config);
  }

  // Combine results from all threads
  std::unordered_map<std::string, uint64_t> total_counts;
  for (const auto &local_map : count_maps) {
    for (const auto &pair : local_map) {
      total_counts[pair.first] += pair.second;
    }
  }
  // convert to DataFrame
  CharacterVector adapters_out(total_counts.size());
  NumericVector counts_out(total_counts.size());
  int i = 0;
  for (const auto &pair : total_counts) {
    adapters_out[i] = pair.first;
    counts_out[i] = static_cast<double>(pair.second);
    i++;
  }
  return DataFrame::create(
      _["adapter"] = adapters_out,
      _["count"]   = counts_out
  );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List rcpp_fraq_summary(const std::vector<std::string> &input,
                             const bool phred33 = true,
                             const int min_overlap = 12,
                             const double max_mismatch_rate = 0.10,
                             const size_t limit = 0,
                             const int nthreads = 1) {
    if (input.size() == 0 || input.size() > 2) {
        throw std::runtime_error("input must contain 1 or 2 fastq files");
    }

    using namespace Rcpp;

    struct Accum {
        std::vector<uint64_t> qsum;
        std::vector<uint64_t> qcnt;
        std::vector<std::unordered_map<char, uint64_t>> base;
        std::unordered_map<int, uint64_t> length_hist;
        std::unordered_map<int, uint64_t> avg_quality_hist;
        uint64_t total_reads = 0;
        uint64_t total_bases = 0;
        uint64_t total_gc = 0;

        inline void ensure_size(size_t need) {
            if (qsum.size() < need) {
                size_t old = qsum.size();
                qsum.resize(need, 0);
                qcnt.resize(need, 0);
                base.resize(need);
                for (size_t i = old; i < need; ++i) base[i].reserve(8);
            }
        }

        inline void process_read(const fraq::Read &r, const bool phred33_flag) {
            const std::string &seq = r.seq;
            const std::string &qual = r.qual;
            const size_t L = seq.size();
            if (L == 0 || qual.size() != L) return;

            total_reads += 1;
            total_bases += L;

            size_t gc_here = 0;
            uint64_t read_q_sum = 0;
            ensure_size(L);

            for (size_t i = 0; i < L; ++i) {
                char b = seq[i];
                if (b >= 'a' && b <= 'z') b = static_cast<char>(b - 'a' + 'A');
                if (b == 'G' || b == 'C') ++gc_here;
                base[i][b] += 1;

                unsigned char qc = static_cast<unsigned char>(qual[i]);
                int q = phred33_flag ? (static_cast<int>(qc) - 33) : (static_cast<int>(qc) - 64);
                if (q < 0) q = 0;
                qsum[i] += static_cast<uint64_t>(q);
                qcnt[i] += 1;
                read_q_sum += static_cast<uint64_t>(q);
            }
            total_gc += gc_here;

            length_hist[static_cast<int>(L)] += 1;

            double avg_quality = static_cast<double>(read_q_sum) / static_cast<double>(L);
            int avg_bin = static_cast<int>(std::llround(avg_quality));
            avg_quality_hist[avg_bin] += 1;
        }
    };

    tbb::enumerable_thread_specific<Accum> acc1_tls;
    tbb::enumerable_thread_specific<Accum> acc2_tls;
    tbb::enumerable_thread_specific<std::unordered_map<int, uint64_t>> insert_tls;

    auto summary_task_single = [&](input_t reads, size_t /*read_index*/) -> fraq::output_t {
        if (reads.size() != 1) {
            throw std::runtime_error("rcpp_fraq_summary(single-end): expected batches of size 1");
        }
        acc1_tls.local().process_read(reads[0], phred33);
        return {};
    };

    auto summary_task_paired = [&](input_t reads, size_t /*read_index*/) -> fraq::output_t {
        if (reads.size() != 2) {
            throw std::runtime_error("rcpp_fraq_summary(paired-end): expected batches of size 2");
        }
        acc1_tls.local().process_read(reads[0], phred33);
        acc2_tls.local().process_read(reads[1], phred33);

        const std::string &s1 = reads[0].seq;
        const std::string &s2 = reads[1].seq;
        if (!s1.empty() && !s2.empty()) {
            std::string s2rc = fraq::reverse_complement(s2);
            int dummy_mismatches = 0;
            size_t ov = fraq_calculate_overlap(s1, s2rc, min_overlap, max_mismatch_rate, dummy_mismatches, false);
            if (ov > 0) {
                int ins = static_cast<int>(s1.size() + s2.size() - ov);
                insert_tls.local()[ins] += 1;
            }
        }
        return {};
    };

    if (input.size() == 1) {
        auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false, limit);
        fraq_run(input, summary_task_single, nthreads, config);
    } else {
        auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false, limit);
        fraq_run(input, summary_task_paired, nthreads, config);
    }

    auto merge_accum = [](const tbb::enumerable_thread_specific<Accum> &tls) -> Accum {
        Accum out;
        for (const auto &A : tls) {
            if (out.qsum.size() < A.qsum.size()) {
                size_t old = out.qsum.size();
                out.qsum.resize(A.qsum.size(), 0);
                out.qcnt.resize(A.qcnt.size(), 0);
                out.base.resize(A.base.size());
                for (size_t i = old; i < A.base.size(); ++i) out.base[i].reserve(8);
            }
            for (size_t i = 0; i < A.qsum.size(); ++i) {
                out.qsum[i] += A.qsum[i];
                out.qcnt[i] += A.qcnt[i];
                auto &dst = out.base[i];
                const auto &src = A.base[i];
                for (const auto &kv : src) dst[kv.first] += kv.second;
            }
            for (const auto &kv : A.length_hist) out.length_hist[kv.first] += kv.second;
            out.total_reads += A.total_reads;
            out.total_bases += A.total_bases;
            out.total_gc += A.total_gc;
            for (const auto &kv : A.avg_quality_hist) out.avg_quality_hist[kv.first] += kv.second;
        }
        return out;
    };

    Accum A1 = merge_accum(acc1_tls);
    Accum A2;
    if (input.size() == 2) A2 = merge_accum(acc2_tls);

    std::unordered_map<int, uint64_t> insert_hist;
    if (input.size() == 2) {
        for (const auto &im : insert_tls) {
            for (const auto &kv : im) insert_hist[kv.first] += kv.second;
        }
    }

    auto make_basic_df = [](const Accum &A) -> Rcpp::DataFrame {
        double gc_pct = A.total_bases ? (100.0 * static_cast<double>(A.total_gc) / static_cast<double>(A.total_bases)) : NA_REAL;
        int min_len = 0;
        int max_len = 0;
        uint64_t sum_len = 0;
        uint64_t sum_cnt = 0;
        if (!A.length_hist.empty()) {
            min_len = A.length_hist.begin()->first;
            max_len = min_len;
            for (const auto &kv : A.length_hist) {
                if (kv.first < min_len) min_len = kv.first;
                if (kv.first > max_len) max_len = kv.first;
                sum_len += static_cast<uint64_t>(kv.first) * kv.second;
                sum_cnt += kv.second;
            }
        }
        double mean_len = sum_cnt ? static_cast<double>(sum_len) / static_cast<double>(sum_cnt) : NA_REAL;
        return Rcpp::DataFrame::create(
            Rcpp::Named("total_sequences") = Rcpp::NumericVector::create(static_cast<double>(A.total_reads)),
            Rcpp::Named("total_bases") = Rcpp::NumericVector::create(static_cast<double>(A.total_bases)),
            Rcpp::Named("seq_len_min") = Rcpp::NumericVector::create(static_cast<double>(min_len)),
            Rcpp::Named("seq_len_mean") = mean_len,
            Rcpp::Named("seq_len_max") = Rcpp::NumericVector::create(static_cast<double>(max_len)),
            Rcpp::Named("gc_percent") = gc_pct
        );
    };

    auto make_quality_df = [](const Accum &A) -> Rcpp::DataFrame {
        const int P = static_cast<int>(A.qsum.size());
        Rcpp::IntegerVector pos(P);
        Rcpp::NumericVector mean_q(P);
        Rcpp::IntegerVector n_q(P);
        for (int i = 0; i < P; ++i) {
            pos[i] = i + 1;
            uint64_t cnt = A.qcnt[static_cast<size_t>(i)];
            n_q[i] = static_cast<int>(cnt);
            mean_q[i] = cnt ? static_cast<double>(A.qsum[static_cast<size_t>(i)]) / static_cast<double>(cnt) : NA_REAL;
        }
        return Rcpp::DataFrame::create(
            Rcpp::Named("position") = pos,
            Rcpp::Named("mean_q") = mean_q,
            Rcpp::Named("n") = n_q
        );
    };

    auto make_base_df_long = [](const Accum &A) -> Rcpp::DataFrame {
        size_t rows = 0;
        for (const auto &m : A.base) rows += m.size();
        Rcpp::IntegerVector position(rows);
        Rcpp::CharacterVector base(rows);
        Rcpp::NumericVector count(rows);
        Rcpp::NumericVector fraction(rows);
        size_t r = 0;
        for (size_t i = 0; i < A.base.size(); ++i) {
            double tot = 0.0;
            for (const auto &kv : A.base[i]) tot += static_cast<double>(kv.second);
            if (tot <= 0.0) tot = 1.0;
            for (const auto &kv : A.base[i]) {
                position[r] = static_cast<int>(i + 1);
                base[r] = std::string(1, kv.first);
                count[r] = static_cast<double>(kv.second);
                fraction[r] = static_cast<double>(kv.second) / tot;
                ++r;
            }
        }
        return Rcpp::DataFrame::create(
            Rcpp::Named("position") = position,
            Rcpp::Named("base") = base,
            Rcpp::Named("count") = count,
            Rcpp::Named("fraction") = fraction
        );
    };

    auto make_length_df = [](const Accum &A) -> Rcpp::DataFrame {
        Rcpp::IntegerVector len(A.length_hist.size());
        Rcpp::NumericVector cnt(A.length_hist.size());
        int i = 0;
        for (const auto &kv : A.length_hist) {
            len[i] = kv.first;
            cnt[i] = static_cast<double>(kv.second);
            ++i;
        }
        return Rcpp::DataFrame::create(
            Rcpp::Named("length") = len,
            Rcpp::Named("count") = cnt
        );
    };

    auto make_avg_quality_df = [](const Accum &A) -> Rcpp::DataFrame {
        Rcpp::IntegerVector avg(A.avg_quality_hist.size());
        Rcpp::NumericVector cnt(A.avg_quality_hist.size());
        int i = 0;
        for (const auto &kv : A.avg_quality_hist) {
            avg[i] = kv.first;
            cnt[i] = static_cast<double>(kv.second);
            ++i;
        }
        return Rcpp::DataFrame::create(
            Rcpp::Named("avg_quality") = avg,
            Rcpp::Named("count") = cnt
        );
    };

    Rcpp::List out;
    out["basic_stats_R1"] = make_basic_df(A1);
    out["per_base_quality_R1"] = make_quality_df(A1);
    out["per_base_content_R1"] = make_base_df_long(A1);
    out["length_distribution_R1"] = make_length_df(A1);
    out["avg_read_quality_R1"] = make_avg_quality_df(A1);

    if (input.size() == 2) {
        out["basic_stats_R2"] = make_basic_df(A2);
        out["per_base_quality_R2"] = make_quality_df(A2);
        out["per_base_content_R2"] = make_base_df_long(A2);
        out["length_distribution_R2"] = make_length_df(A2);
        out["avg_read_quality_R2"] = make_avg_quality_df(A2);

        if (!insert_hist.empty()) {
            Rcpp::IntegerVector isize(insert_hist.size());
            Rcpp::NumericVector icount(insert_hist.size());
            int i = 0;
            for (const auto &kv : insert_hist) {
                isize[i] = kv.first;
                icount[i] = static_cast<double>(kv.second);
                ++i;
            }
            out["insert_size"] = Rcpp::DataFrame::create(
                Rcpp::Named("insert_size") = isize,
                Rcpp::Named("count") = icount
            );
        }
    }

    return out;
}

// legacy entry point without limit (keeps compatibility with existing registrations)
Rcpp::List rcpp_fraq_summary(const std::vector<std::string> &input,
                             const bool phred33,
                             const int min_overlap,
                             const double max_mismatch_rate,
                             const int nthreads) {
    return rcpp_fraq_summary(input, phred33, min_overlap, max_mismatch_rate, static_cast<size_t>(0), nthreads);
}


// [[Rcpp::export(rng=false)]]
DataFrame rcpp_fraq_align(
    CharacterVector query,
    CharacterVector target,
    int           max_distance    = 2147483647, // fraq::MAX_INT
    std::string   ambiguity_base  = "",
    std::string   boundary        = "contains",
    std::string   distance_metric = "lv"
) {
    int n = query.size();
    if (target.size() != n) {
        throw std::runtime_error("`query` and `target` must have the same length");
    }
    if( boundary != "contains" && boundary != "starts" && boundary != "global") {
        throw std::runtime_error("`boundary` must be either 'contains', 'starts' or 'global'");
    }
    if (distance_metric != "lv" && distance_metric != "levenshtein" && distance_metric != "hm" && distance_metric != "hamming") {
        throw std::runtime_error("`distance_metric` must be either 'lv' or 'hm'");
    }
    char ab = ambiguity_base.empty() ? '\0' : ambiguity_base[0];

    IntegerVector starts(n), ends(n), dists(n);
    for (int i = 0; i < n; ++i) {
      std::string q = as<std::string>(query[i]);
      std::string t = as<std::string>(target[i]);
      fraq::AlignResult res;

      if (distance_metric == "lv" || distance_metric == "levenshtein") {
        if (boundary == "starts") {
            res = fraq_lv_starts(q, t, max_distance, ab);
        } else if(boundary == "contains") {
            res = fraq_lv_contains(q, t, max_distance, ab);
        } else {
            res = fraq_lv_global(q, t, max_distance, ab);
        }
      } else {
        if (boundary == "starts") {
            res = fraq_hm_starts(q, t, max_distance, ab);
        } else if(boundary == "contains") {
            res = fraq_hm_contains(q, t, max_distance, ab);
        } else {
            res = fraq_hm_global(q, t, max_distance, ab);
        }
      }

      starts[i] = static_cast<int>(res.start);
      ends[i]   = static_cast<int>(res.end);
      dists[i]  = res.distance;
    }

    return DataFrame::create(
        _["start"]    = starts,
        _["end"]      = ends,
        _["distance"] = dists
    );
}

// kernel intended for internal testing
// [[Rcpp::export(rng=false)]]
void rcpp_fraq_wait(std::vector<std::string> input,
                    const std::vector<std::string> &output,
                    const int sleep_ms,
                    const int nthreads = 1) {
  if (input.size() == 0 || input.size() != output.size()) {
    throw std::runtime_error("input and output must be character vectors of the same length > 0");
  }
  if (sleep_ms < 0) {
    throw std::runtime_error("sleep_ms must be >= 0");
  }
  auto wait_task = [&](input_t reads, size_t index) -> output_t {
    (void)index;
    if (sleep_ms > 0) {
      std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    }
    return zip(output, std::move(reads));
  };
  auto config = fraq::FraqRunConfig(BLOCKSIZE, FRAQ_COMPRESS_LEVEL, ZSTD_COMPRESS_LEVEL, GZIP_COMPRESS_LEVEL, false);
  fraq_run(input, wait_task, nthreads, config);
}


// [[Rcpp::export(rng=false)]]
DataFrame rcpp_fraq_mem_list() {
    std::vector<std::string> keys;
    std::vector<double> counts;
    keys.reserve(mem_store.size());
    counts.reserve(mem_store.size());
    for (const auto &kv : mem_store) {
        size_t total_reads = 0;
        for (const auto &blk : kv.second) {
            total_reads += static_cast<size_t>(blk.num_reads);
        }
        keys.push_back(kv.first);
        counts.push_back(static_cast<double>(total_reads));
    }
    CharacterVector key_vec(keys.begin(), keys.end());
    NumericVector count_vec(counts.begin(), counts.end());
    return DataFrame::create(
        _["mem_key"] = key_vec,
        _["n_reads"] = count_vec,
        _["stringsAsFactors"] = false
    );
}

// [[Rcpp::export(rng=false)]]
bool rcpp_fraq_mem_remove(const std::string &mem_key) {
    auto it = mem_store.find(mem_key);
    if (it == mem_store.end()) {
        return false;
    }
    mem_store.erase(it);
    return true;
}

// [[Rcpp::export(rng=false)]]
bool fraq_fifo_supported() {
#if FRAQ_HAS_FIFO
    return true;
#else
    return false;
#endif
}
