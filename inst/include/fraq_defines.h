// header is used both internally and externally
// header contains all user facing types, definitions and simple inline functions

#ifndef __FRAQ_DEFINES_H
#define __FRAQ_DEFINES_H

#include <functional>
#include <string>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
#include <stdexcept>
#include <string_view>
#include <sys/stat.h>

#if !defined(S_IFIFO) && defined(_S_IFIFO)
#  define S_IFIFO _S_IFIFO
#endif

#if !defined(FRAQ_HAS_FIFO)
#  if (defined(__unix__) || defined(__APPLE__)) && defined(S_IFIFO)
#    define FRAQ_HAS_FIFO 1
#  else
#    define FRAQ_HAS_FIFO 0
#  endif
#endif
	

namespace fraq {

struct Read {
  std::string name;
  std::string seq;
  std::string qual;
};

inline void validate_read_storage_invariants(const Read &read) {
    if (read.seq.size() != read.qual.size()) {
        throw std::runtime_error("read sequence and quality lengths differ");
    }
}

// tuple implementation
template<size_t I> auto& get(Read& r) noexcept {
    if constexpr (I == 0) return r.name;
    if constexpr (I == 1) return r.seq;
    return r.qual;
}
template<size_t I> const auto& get(const Read& r) noexcept {
    if constexpr (I == 0) return r.name;
    if constexpr (I == 1) return r.seq;
    return r.qual;
}
template<size_t I> auto&& get(Read&& r) noexcept {
    if constexpr (I == 0) return std::move(r.name);
    if constexpr (I == 1) return std::move(r.seq);
    return std::move(r.qual);
}

struct AlignResult {
  size_t start;
  size_t end;
  int distance;
};

using process_task_t = std::function< std::vector< std::pair<std::string, Read> >(std::vector<Read>, std::size_t) >;
inline static const std::string NO_OUTPUT = ""; // no constexpr for std::string until C++20
constexpr int MAX_INT = std::numeric_limits<int>::max();
using input_t = std::vector<Read>;
using outpair_t = std::pair<std::string, Read>;
using output_t = std::vector<outpair_t>;

inline std::string reverse_complement(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.size());

    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        char c;
        switch (*it) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            case 'a': c = 't'; break;
            case 'c': c = 'g'; break;
            case 'g': c = 'c'; break;
            case 't': c = 'a'; break;
            default:  c = *it; break;
        }
        rc.push_back(c);
    }

    return rc;
}

inline void trim_read(Read &read, const size_t start, const size_t end) {
    const size_t seqlen = read.seq.size();
    if (end > seqlen) {
        throw std::runtime_error("trim_read end past sequence length");
    }
    if (start > end) {
        throw std::runtime_error("trim_read start greater than end");
    }
    const size_t len = end - start;
    read.seq  = read.seq.substr(start, len);
    read.qual = read.qual.substr(start, len);
}

template<class A, class B>
inline std::vector<std::pair<A,B>> zip(std::vector<A> va, std::vector<B> vb) {
    size_t n = std::min(va.size(), vb.size());
    std::vector<std::pair<A,B>> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i)
        out.emplace_back(std::move(va[i]), std::move(vb[i]));
    return out;
}

inline bool ends_with(const std::string &s, const std::string &suffix) {
    return s.size() >= suffix.size()
        && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

struct FraqRunConfig {
    // size_t blocksize
    // int compress values for fraq/mem, zstd, gzip
    // default parameters and validate values are within acceptable ranges
    size_t blocksize = 65535L;
    int fraq_compress_level = 3;
    int zstd_compress_level = 3;
    int gzip_compress_level = 6;
    bool serial_kernel = false; // use serial kernel or unlimited kernel
    size_t limit = 0; // 0 means no limit
    bool (*interrupt_fn)(void* ctx) = nullptr;
    void* interrupt_ctx = nullptr;
    FraqRunConfig() {}
    FraqRunConfig(size_t blocksize, int fraq_compress_level, int zstd_compress_level, int gzip_compress_level, bool serial_kernel, size_t limit = 0)
      : blocksize(blocksize),
        fraq_compress_level(fraq_compress_level),
        zstd_compress_level(zstd_compress_level),
        gzip_compress_level(gzip_compress_level),
        serial_kernel(serial_kernel),
        limit(limit)
    {
        if (blocksize == 0 || blocksize > 65535L) {
            throw std::runtime_error("blocksize must be between 1 and 65535");
        }
        if (fraq_compress_level < 1 || fraq_compress_level > 22) {
            throw std::runtime_error("fraq compress level must be between 1 and 22");
        }
        if (zstd_compress_level < 1 || zstd_compress_level > 22) {
            throw std::runtime_error("zstd compress level must be between 1 and 22");
        }
        if (gzip_compress_level < 1 || gzip_compress_level > 9) {
            throw std::runtime_error("gzip compress level must be between 1 and 9");
        }
    }
};

}; // namespace fraq


// tuple interface for fraq::Read
namespace std {
    template<> struct tuple_size<fraq::Read> : integral_constant<size_t,3> {};
    template<> struct tuple_element<0, fraq::Read> { using type = std::string; };
    template<> struct tuple_element<1, fraq::Read> { using type = std::string; };
    template<> struct tuple_element<2, fraq::Read> { using type = std::string; };
}


#endif // include guard
