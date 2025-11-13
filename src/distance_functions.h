#ifndef __FRAQ_DISTANCE_FUNCTIONS_H
#define __FRAQ_DISTANCE_FUNCTIONS_H

#include <string>
#include <string_view>
#include <array>
#include <vector>
#include "fraq_defines.h"
#include "edlib_external.h"

int fraq_global_hamming_distance(const char* a, const char* b, const size_t n, const int max_distance, const char ambiguity_base) noexcept;
fraq::AlignResult fraq_hm_starts(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') noexcept;
fraq::AlignResult fraq_hm_contains(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') noexcept;
fraq::AlignResult fraq_hm_global(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0') noexcept;
fraq::AlignResult fraq_lv_starts(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0');
fraq::AlignResult fraq_lv_contains(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0');
fraq::AlignResult fraq_lv_global(std::string_view query, std::string_view target, int max_distance = fraq::MAX_INT, char ambiguity_base = '\0');

#if defined(__has_builtin)
  #if __has_builtin(__builtin_popcount)
    #define HAVE_BUILTIN_POPCOUNT 1
  #else
    #define HAVE_BUILTIN_POPCOUNT 0
  #endif
#else
  #if defined(__GNUC__) || defined(__clang__)
    #define HAVE_BUILTIN_POPCOUNT 1
  #else
    #define HAVE_BUILTIN_POPCOUNT 0
  #endif
#endif

#if HAVE_BUILTIN_POPCOUNT && defined(__SSE2__)
    #define USE_SIMD_HAMMING
#endif

#ifdef USE_SIMD_HAMMING
  #include <emmintrin.h>
#endif

int fraq_global_hamming_distance(
    const char* a,
    const char* b,
    const size_t n,
    const int max_distance,
    const char ambiguity_base) noexcept {
    int dist = 0;

  #ifdef USE_SIMD_HAMMING
    const __m128i all_ones = _mm_set1_epi8(static_cast<char>(-1));
    const __m128i amb_vec  = _mm_set1_epi8(ambiguity_base);

    size_t i = 0;
    size_t simd_end = (n / 16) * 16;
    for (; i < simd_end; i += 16) {
        __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a + i));
        __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b + i));

        __m128i eq      = _mm_cmpeq_epi8(va, vb);
        __m128i m_qamb  = _mm_cmpeq_epi8(va, amb_vec);
        __m128i m_tamb  = _mm_cmpeq_epi8(vb, amb_vec);
        __m128i match   = _mm_or_si128(eq, _mm_or_si128(m_qamb, m_tamb));

        __m128i mismatch = _mm_andnot_si128(match, all_ones);
        int     mask     = _mm_movemask_epi8(mismatch);

        dist += __builtin_popcount(static_cast<unsigned>(mask));
        if (dist > max_distance) return fraq::MAX_INT;
    }

    for (; i < n; ++i) {
        char qc = a[i], tc = b[i];
        if (qc != tc && qc != ambiguity_base && tc != ambiguity_base) {
            if (++dist > max_distance) return fraq::MAX_INT;
        }
    }
    return dist;

  #else
    for (size_t i = 0; i < n; ++i) {
        char qc = a[i], tc = b[i];
        if (qc != tc && qc != ambiguity_base && tc != ambiguity_base) {
            ++dist;
        }
    }
    return (dist <= max_distance ? dist : fraq::MAX_INT);
  #endif
}

fraq::AlignResult fraq_hm_starts(std::string_view query, std::string_view target,
                      int max_distance, char ambiguity_base) noexcept {
    const size_t qlen = query.size();
    const size_t tlen = target.size();
    if (qlen > tlen) {
        return {0, 0, fraq::MAX_INT};
    }
    if (qlen == 0) {
        return {0, 0, 0};
    }

    fraq::AlignResult out = {0, 0, fraq::MAX_INT};
    int result = fraq_global_hamming_distance(
        query.data(), target.data(), qlen,
        max_distance, ambiguity_base
    );
    if (result <= max_distance) {
        out.end = qlen - 1;
        out.distance = result;
    }
    return out;
}

fraq::AlignResult fraq_hm_contains(
    std::string_view query,
    std::string_view target,
    int                max_distance,
    char               ambiguity_base
) noexcept {
    const size_t qlen = query.size();
    const size_t tlen = target.size();
    if (qlen > tlen) {
        return {0, 0, fraq::MAX_INT};
    }
    if (qlen == 0) {
        return {0, 0, 0};
    }

    const char* q = query.data();
    const char* t = target.data();

    size_t best_i = 0;
    int best_distance = fraq::MAX_INT;
    for (size_t i = 0; i + qlen <= tlen; ++i) {
        int d = fraq_global_hamming_distance(q, t + i, qlen, max_distance, ambiguity_base);
        if(d == 0) {
            return {i, i + qlen - 1, 0};
        }
        if(d < best_distance) {
            best_distance = d;
            best_i = i;
        }
    }
    fraq::AlignResult out = {0, 0, fraq::MAX_INT};
    if (best_distance <= max_distance) {
        out.start    = best_i;
        out.end      = best_i + qlen - 1;
        out.distance = best_distance;
        return out;
    }
    return out;
}

fraq::AlignResult fraq_hm_global(std::string_view query, std::string_view target,
                                 int max_distance, char ambiguity_base) noexcept {
    const size_t qlen = query.size();
    const size_t tlen = target.size();
    fraq::AlignResult out = {0, 0, fraq::MAX_INT};

    if (qlen != tlen) return out;
    if (qlen == 0) { out.distance = 0; return out; }

    int d = fraq_global_hamming_distance(query.data(), target.data(), qlen,
                                         max_distance, ambiguity_base);
    if (d <= max_distance) {
        out.start = 0;
        out.end = qlen - 1;
        out.distance = d;
    }
    return out;
}

inline std::vector<EdlibEqualityPair> build_equality_pairs(
    std::string_view query,
    std::string_view target,
    char ambiguity_base
) {
  if (ambiguity_base == '\0') {
    return {};
  }

  const unsigned char amb = static_cast<unsigned char>(ambiguity_base);
  std::array<bool, 256> seen{};

  for (char c : query) {
    unsigned char uc = static_cast<unsigned char>(c);
    if (uc != amb) seen[uc] = true;
  }
  for (char c : target) {
    unsigned char uc = static_cast<unsigned char>(c);
    if (uc != amb) seen[uc] = true;
  }

  size_t count = 0;
  for (bool flag : seen) {
    if (flag) ++count;
  }

  std::vector<EdlibEqualityPair> eqs;
  eqs.reserve(count);

  for (int i = 0; i < 256; ++i) {
    if (seen[i]) {
      eqs.emplace_back(EdlibEqualityPair{ ambiguity_base, static_cast<char>(i) });
    }
  }

  return eqs;
}

fraq::AlignResult fraq_lv_starts(std::string_view query, std::string_view target,
                      int max_distance, char ambiguity_base) {
  auto eqs = build_equality_pairs(query, target, ambiguity_base);
  int eqs_len = static_cast<int>(eqs.size());
  const EdlibEqualityPair* eqs_ptr = eqs_len ? eqs.data() : nullptr;

  EdlibAlignConfig config = edlib_new_align_config(
    (max_distance == fraq::MAX_INT ? -1 : max_distance),
    EDLIB_MODE_SHW,
    EDLIB_TASK_LOC,
    eqs_ptr,
    eqs_len
  );

  EdlibAlignResult result = edlib_align(
    query.data(),  static_cast<int>(query.size()),
    target.data(), static_cast<int>(target.size()),
    config
  );
  if(result.status) throw(std::runtime_error("edlib_align returned an error"));

  fraq::AlignResult out{0, 0, fraq::MAX_INT};
  if(result.editDistance == -1 || result.numLocations == 0 || result.endLocations[0] == -1) return out;
  out.start = result.startLocations[0];
  out.end   = result.endLocations[0];
  out.distance = result.editDistance;

  edlib_free_align_result(result);

  return out;
}

fraq::AlignResult fraq_lv_contains(std::string_view query, std::string_view target,
                      int max_distance, char ambiguity_base) {
  auto eqs = build_equality_pairs(query, target, ambiguity_base);
  int eqs_len = static_cast<int>(eqs.size());
  const EdlibEqualityPair* eqs_ptr = eqs_len ? eqs.data() : nullptr;

  EdlibAlignConfig config = edlib_new_align_config(
    (max_distance == fraq::MAX_INT ? -1 : max_distance),
    EDLIB_MODE_HW,
    EDLIB_TASK_LOC,
    eqs_ptr,
    eqs_len
  );

  EdlibAlignResult result = edlib_align(
    query.data(),  static_cast<int>(query.size()),
    target.data(), static_cast<int>(target.size()),
    config
  );
  if(result.status) throw(std::runtime_error("edlib_align returned an error"));

  fraq::AlignResult out{0, 0, fraq::MAX_INT};
  if(result.editDistance == -1 || result.numLocations == 0 || result.endLocations[0] == -1) return out;
  out.start = result.startLocations[0];
  out.end   = result.endLocations[0];
  out.distance = result.editDistance;

  edlib_free_align_result(result);

  return out;
}

fraq::AlignResult fraq_lv_global(std::string_view query, std::string_view target,
                           int max_distance, char ambiguity_base) {
    auto eqs = build_equality_pairs(query, target, ambiguity_base);
    int eqs_len = static_cast<int>(eqs.size());
    const EdlibEqualityPair* eqs_ptr = eqs_len ? eqs.data() : nullptr;

    EdlibAlignConfig config = edlib_new_align_config(
        (max_distance == fraq::MAX_INT ? -1 : max_distance),
        EDLIB_MODE_NW,
        EDLIB_TASK_LOC,
        eqs_ptr,
        eqs_len
    );

    EdlibAlignResult result = edlib_align(
        query.data(), static_cast<int>(query.size()),
        target.data(), static_cast<int>(target.size()),
        config
    );
    if (result.status) throw(std::runtime_error("edlib_align returned an error"));

    fraq::AlignResult out{0, 0, fraq::MAX_INT};
    if (result.editDistance != -1 && result.numLocations > 0 && result.endLocations[0] != -1) {
        out.start = result.startLocations[0];
        out.end = result.endLocations[0];
        out.distance = result.editDistance;
    }

    edlib_free_align_result(result);
    return out;
}

#endif // include guard
