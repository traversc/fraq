#ifndef _FRAQ_CONCAT_SERIAL_H_
#define _FRAQ_CONCAT_SERIAL_H_

#include "fraq_run_serial.h"

namespace fraq_internal {

inline void serial_fraq_concat(
    const std::vector<std::string> &input_files,
    const std::string &output_file,
    const FraqRunConfig &config,
    std::unordered_map<std::string, std::vector<CompressedReadBlock>> &mem_store) {
  if (input_files.empty()) {
    throw std::runtime_error("fraq_concat requires at least one input file");
  }
  if (output_file.empty()) {
    throw std::runtime_error("fraq_concat output path must be non-empty");
  }

  std::vector<MultiReader> readers;
  readers.reserve(input_files.size());
  for (const auto &file : input_files) {
    readers.emplace_back(make_serial_reader(file, mem_store));
  }

  SerialOutputRouter router(config, mem_store);
  try {
    for (auto &reader : readers) {
      while (true) {
        std::vector<Read> block;
        reader.read_block(block, config.blocksize);
        if (block.empty()) {
          break;
        }
        auto processed = std::make_shared<ProcessedBlock>();
        processed->index = 0;
        processed->bins[output_file] = std::move(block);
        router.write_processed(processed);
      }
    }
    router.finish();
  } catch (...) {
    router.cleanup_outputs();
    throw;
  }
}

} // namespace fraq_internal

#endif // _FRAQ_CONCAT_SERIAL_H_
