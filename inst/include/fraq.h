// this header is used ONLY externally

// Notes on code organization and namespaces
// namespace fraq: external facing, everything needed for user to implement new fraq kernel
// namespace fraq_internal: internal implementation details, IO, multitheading
// src/fraq_classes.h fraq_internal classes
// src/fraq_functions.cpp
//    Contains fraq_run function
//    Contains pre-compiled kernels
//    Contains helper functions (hamming, levenshtein distance)
// fraq_defines.h all user facing types and definitions
// fraq.h
//    Exports fraq_run (becomes fraq::run for external consistency) and helper functions
//    Exports helper functions
//    extern "C" is not necessary. But we do need inline because RcppExports.cpp auto-imports fraq.h so we could run into ODR issues without it 

#ifndef __FRAQ_H
#define __FRAQ_H

#include "fraq_defines.h"
#include <R_ext/Rdynload.h>

namespace fraq {
  inline void run(const std::vector<std::string> &input_files, process_task_t task, const int nthreads = 1, FraqRunConfig config = FraqRunConfig{}) {
    static void (*fun)(const std::vector<std::string>&, process_task_t, const int, FraqRunConfig) = 
      (void (*)(const std::vector<std::string>&, process_task_t, const int, FraqRunConfig)) R_GetCCallable("fraq", "fraq_run");
    fun(input_files, task, nthreads, config);
  }

  inline void concat(const std::vector<std::string> &input_files,
                     const std::string &output_file,
                     const int nthreads = 1,
                     FraqRunConfig config = FraqRunConfig{}) {
    static void (*fun)(const std::vector<std::string>&, const std::string&, const int, const FraqRunConfig) =
      (void (*)(const std::vector<std::string>&, const std::string&, const int, const FraqRunConfig)) R_GetCCallable("fraq", "fraq_concat");
    fun(input_files, output_file, nthreads, config);
  }
} // namespace fraq

namespace fraq {
  using align_fn_t = AlignResult (*)(std::string_view, std::string_view, int, char);

  inline AlignResult hm_starts(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') noexcept {
      static align_fn_t fun = reinterpret_cast<align_fn_t>(R_GetCCallable("fraq", "fraq_hm_starts"));
      return fun(query, target, max_distance, ambiguity_base);
  }

  inline AlignResult hm_contains(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') noexcept {
      static align_fn_t fun = reinterpret_cast<align_fn_t>(R_GetCCallable("fraq", "fraq_hm_contains"));
      return fun(query, target, max_distance, ambiguity_base);
  }

  inline AlignResult hm_global(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') noexcept {
      static align_fn_t fun = reinterpret_cast<align_fn_t>(R_GetCCallable("fraq", "fraq_hm_global"));
      return fun(query, target, max_distance, ambiguity_base);
  }

  inline AlignResult lv_starts(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') {
      static align_fn_t fun = reinterpret_cast<align_fn_t>(R_GetCCallable("fraq", "fraq_lv_starts"));
      return fun(query, target, max_distance, ambiguity_base);
  }

  inline AlignResult lv_contains(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') {
      static align_fn_t fun = reinterpret_cast<align_fn_t>(R_GetCCallable("fraq", "fraq_lv_contains"));
      return fun(query, target, max_distance, ambiguity_base);
  }

  inline AlignResult lv_global(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') {
      static align_fn_t fun = reinterpret_cast<align_fn_t>(R_GetCCallable("fraq", "fraq_lv_global"));
      return fun(query, target, max_distance, ambiguity_base);
  }
} // namespace fraq


#endif // include guard
