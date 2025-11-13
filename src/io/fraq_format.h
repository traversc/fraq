#ifndef __FRAQ_FORMAT_H
#define __FRAQ_FORMAT_H

#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <limits>
#include <fstream>
#include <iostream>
#include <zstd.h>
#include "fraq_defines.h"

namespace fraq_internal {

using fraq::Read;

constexpr size_t MAX_BLOCKSIZE = 65535; // max 16-bit block size

// standard nucleotides for bit packing
static const auto fraq_standard_bases = [](){
  std::array<uint8_t,256> m;
  m.fill(0xFF); // 0xFF means invalid
  m['A'] = 0; m['C'] = 1; m['G'] = 2; m['T'] = 3;
  m['R'] = 4; m['Y'] = 5; m['S'] = 6; m['W'] = 7;
  m['K'] = 8; m['M'] = 9; m['B'] = 10; m['D'] = 11;
  m['H'] = 12; m['V'] = 13; m['N'] = 14; m['U'] = 15;
  return m;
}();

// max compress level to do nucleotide bit packing
// This is a heuristic, at max compress_level NOT packing results in smaller output
constexpr int MAX_COMPRESS_LEVEL_FOR_PACKING = 19;

// dynamic array that holds uint8_t, uint16_t, uint32_t or uint64_t
// We use memcpy to set and get values, so we dont have to worry about alignment. 
struct DynArray {
  enum class ElemSize : uint8_t { B8 = 0, B16 = 1, B32 = 2, B64 = 3 };
  
  ElemSize          sz;
  std::vector<char> buf;
  
  static constexpr std::size_t bytes_per_element(ElemSize s) {
    switch (s) {
    case ElemSize::B8:  return 1;
    case ElemSize::B16: return 2;
    case ElemSize::B32: return 4;
    case ElemSize::B64: return 8;
    }
    return 1;
  }
  
  static ElemSize choose_size(std::size_t max_len) {
    if (max_len <= std::numeric_limits<uint8_t >::max()) return ElemSize::B8;
    if (max_len <= std::numeric_limits<uint16_t>::max()) return ElemSize::B16;
    if (max_len <= std::numeric_limits<uint32_t>::max()) return ElemSize::B32;
    return ElemSize::B64;
  }
  
  static DynArray make(std::size_t n, std::size_t max_len) {
    return DynArray(n, choose_size(max_len));
  }
  
  DynArray(std::size_t n, ElemSize s) : sz(s), buf(n * bytes_per_element(s), 0) {}
  
  uint64_t get(std::size_t i) const {
    const std::size_t stride = bytes_per_element(sz);
    uint64_t v = 0;
    std::memcpy(&v, buf.data() + i * stride, stride);
    return v;
  }
  
  void set(std::size_t i, uint64_t v) {
    const std::size_t stride = bytes_per_element(sz);
    std::memcpy(buf.data() + i * stride, &v, stride);
  }
  
  struct Proxy {
    DynArray&   arr;
    std::size_t idx;
    Proxy(DynArray& a, std::size_t i) : arr(a), idx(i) {}
    Proxy& operator=(uint64_t v) { arr.set(idx, v); return *this; }
    operator uint64_t() const { return arr.get(idx); }
  };
  
  Proxy    operator[](std::size_t i)       { return Proxy(*this, i); }
  uint64_t operator[](std::size_t i) const { return get(i); }
};

bool is_big_endian() {
  uint32_t x = 0x01020304u;
  unsigned char b[4];
  std::memcpy(b, &x, 4);
  return b[0] == 1;
}

std::vector<char> do_zstd_compress(const std::vector<char>& buf, const int compress_level){
  size_t bound = ZSTD_compressBound(buf.size());
  std::vector<char> out(bound);
  size_t csize = ZSTD_compress(out.data(), bound, buf.data(), buf.size(), compress_level);
  if (ZSTD_isError(csize)) {
    throw std::runtime_error(ZSTD_getErrorName(csize));
  }
  out.resize(csize);
  return out;
};


// Container for compressed block data.
struct CompressedReadBlock {
  std::string         name_prefix;
  std::vector<char>   compressed_name_lengths;
  std::vector<char>   compressed_names;
  std::vector<char>   compressed_seq_lengths;
  std::vector<char>   compressed_seqs;
  std::vector<char>   compressed_quals;
  uint64_t            uncompressed_names_size; // uncompressed size
  uint64_t            uncompressed_seqs_size; // uncompressed size
  uint16_t            num_reads;
  DynArray::ElemSize  name_len_elem_size;
  DynArray::ElemSize  seq_len_elem_size;
  bool                use_bit_pack;
};


////////////////////////////////////////////////////////////////////////////////
// write

CompressedReadBlock compress_reads(const Read * reads, const size_t n, const int compress_level) {
  if (n > MAX_BLOCKSIZE)
    throw std::runtime_error("block size > MAX_BLOCKSIZE");
  
  // 1) scan for sizes, non-standard bases, and length overflow
  size_t total_name_bytes = 0;
  size_t total_seq_bytes  = 0;
  size_t max_seq_len      = 0;
  size_t max_name_len     = 0;
  
  std::string common_name_prefix = reads[0].name;
  size_t prefix_limit = common_name_prefix.size();
  bool use_bit_pack = compress_level <= MAX_COMPRESS_LEVEL_FOR_PACKING;
  for (size_t i = 0; i < n; ++i) {
    auto & m = reads[i];
    size_t name_len = m.name.size();
    size_t seq_len  = m.seq.size();
    total_name_bytes += name_len;
    total_seq_bytes  += seq_len;
    max_seq_len = std::max(max_seq_len, seq_len);
    max_name_len = std::max(max_name_len, name_len);
    
    prefix_limit = std::min(prefix_limit, name_len);
    auto pr = std::mismatch(common_name_prefix.begin(), common_name_prefix.begin() + prefix_limit, m.name.begin());
    prefix_limit = static_cast<size_t>(std::distance(common_name_prefix.begin(), pr.first));
    
    if (use_bit_pack) {
      for (const unsigned char c : m.seq) {
        if (fraq_standard_bases[c] == 0xFF) { use_bit_pack = false; break; }
      }
    }
  }
  common_name_prefix.resize(prefix_limit);
  

  DynArray name_sizes = DynArray::make(n, max_name_len - common_name_prefix.size());
  DynArray seq_sizes = DynArray::make(n, max_seq_len);
  std::vector<char> names(total_name_bytes - n * common_name_prefix.size() );
  std::vector<char> seqs(use_bit_pack ? ((total_seq_bytes + 1) / 2) : total_seq_bytes);
  std::vector<char> quals(total_seq_bytes);
  
  // 1) write names, seqs (packed or raw), quals
  char *pname = names.data();
  char *pseq  = seqs.data();
  char *pqual = quals.data();
  
  size_t pos = 0; // global base index for packing
  for (size_t i = 0; i < n; ++i) {
    auto const &name = reads[i].name;
    auto const &seq  = reads[i].seq;
    auto const &qual = reads[i].qual;
    
    const size_t name_skip = common_name_prefix.size();
    const size_t name_len  = name.size() - name_skip;
    
    name_sizes.set(i, name_len);
    std::memcpy(pname, name.data() + name_skip, name_len);
    pname += name_len;
    
    seq_sizes.set(i, seq.size());
    if (use_bit_pack) {
      for (char b : seq) {
        const uint8_t code = fraq_standard_bases[static_cast<unsigned char>(b)];
        const size_t j = pos >> 1;
        if ((pos & 1) == 0) {
          pseq[j] = static_cast<char>(code << 4);
        } else {
          pseq[j] = static_cast<char>(static_cast<uint8_t>(pseq[j]) | code);
        }
        ++pos;
      }
    } else {
      std::memcpy(pseq, seq.data(), seq.size());
      pseq += seq.size();
    }
    std::memcpy(pqual, qual.data(), qual.size());
    pqual += qual.size();
  }
  
  // return block
  auto out = CompressedReadBlock{
    std::move(common_name_prefix),
    do_zstd_compress(name_sizes.buf, compress_level),
    do_zstd_compress(names, compress_level),
    do_zstd_compress(seq_sizes.buf, compress_level),
    do_zstd_compress(seqs, compress_level),
    do_zstd_compress(quals, compress_level),
    static_cast<uint64_t>(names.size()),
    static_cast<uint64_t>(seqs.size()),
    static_cast<uint16_t>(n),
    name_sizes.sz,
    seq_sizes.sz,
    use_bit_pack
  };
  return out;
}

CompressedReadBlock compress_reads(const std::vector<Read> & reads, const int compress_level) {
  return compress_reads(reads.data(), reads.size(), compress_level);
}

std::vector<char> generate_block_header(const CompressedReadBlock &b) {
  using ES = DynArray::ElemSize;
  
  const uint64_t v_num_reads              = static_cast<uint64_t>(b.num_reads);
  const uint64_t v_uncompressed_names_sz  = b.uncompressed_names_size;
  const uint64_t v_uncompressed_seqs_sz   = b.uncompressed_seqs_size;
  const uint64_t v_name_prefix_sz         = static_cast<uint64_t>(b.name_prefix.size());
  const uint64_t v_cname_lengths_sz       = static_cast<uint64_t>(b.compressed_name_lengths.size());
  const uint64_t v_cnames_sz              = static_cast<uint64_t>(b.compressed_names.size());
  const uint64_t v_cseq_lengths_sz        = static_cast<uint64_t>(b.compressed_seq_lengths.size());
  const uint64_t v_cseqs_sz               = static_cast<uint64_t>(b.compressed_seqs.size());
  const uint64_t v_cquals_sz              = static_cast<uint64_t>(b.compressed_quals.size());
  
  const ES c_num_reads             = DynArray::choose_size(v_num_reads);
  const ES c_uncompressed_names_sz = DynArray::choose_size(v_uncompressed_names_sz);
  const ES c_uncompressed_seqs_sz  = DynArray::choose_size(v_uncompressed_seqs_sz);
  const ES c_name_prefix_sz        = DynArray::choose_size(v_name_prefix_sz);
  const ES c_cname_lengths_sz      = DynArray::choose_size(v_cname_lengths_sz);
  const ES c_cnames_sz             = DynArray::choose_size(v_cnames_sz);
  const ES c_cseq_lengths_sz       = DynArray::choose_size(v_cseq_lengths_sz);
  const ES c_cseqs_sz              = DynArray::choose_size(v_cseqs_sz);
  const ES c_cquals_sz             = DynArray::choose_size(v_cquals_sz);
  
  uint32_t meta = 0;
  auto pack = [&](uint32_t code, uint32_t shift) {
    meta |= (code & 0x3u) << shift; // 2-bit field (supports B8,B16,B32,B64)
  };
  pack(static_cast<uint32_t>(c_num_reads),              0);
  pack(static_cast<uint32_t>(c_uncompressed_names_sz),  2);
  pack(static_cast<uint32_t>(c_uncompressed_seqs_sz),   4);
  pack(static_cast<uint32_t>(c_name_prefix_sz),         6);
  pack(static_cast<uint32_t>(c_cname_lengths_sz),       8);
  pack(static_cast<uint32_t>(c_cnames_sz),             10);
  pack(static_cast<uint32_t>(c_cseq_lengths_sz),       12);
  pack(static_cast<uint32_t>(c_cseqs_sz),              14);
  pack(static_cast<uint32_t>(c_cquals_sz),             16);
  pack(static_cast<uint32_t>(b.name_len_elem_size),    18);
  pack(static_cast<uint32_t>(b.seq_len_elem_size),     20);
  meta |= (static_cast<uint32_t>(b.use_bit_pack) & 0x1u) << 22;
  
  auto byte_count = [](ES es) -> size_t {
    return es == ES::B8  ? 1 :
    es == ES::B16 ? 2 :
    es == ES::B32 ? 4 : 8; // B64 or future
  };
  
  size_t total_bytes = 4;
  total_bytes += byte_count(c_num_reads);
  total_bytes += byte_count(c_uncompressed_names_sz);
  total_bytes += byte_count(c_uncompressed_seqs_sz);
  total_bytes += byte_count(c_name_prefix_sz);
  total_bytes += byte_count(c_cname_lengths_sz);
  total_bytes += byte_count(c_cnames_sz);
  total_bytes += byte_count(c_cseq_lengths_sz);
  total_bytes += byte_count(c_cseqs_sz);
  total_bytes += byte_count(c_cquals_sz);
  
  std::vector<char> out;
  out.reserve(total_bytes);
  
  out.resize(4);
  std::memcpy(out.data(), &meta, 4);
  
  auto append_val = [&](uint64_t v, ES es) {
    const size_t n = byte_count(es);
    const size_t off = out.size();
    out.resize(off + n);
    std::memcpy(out.data() + off, &v, n);
  };
  
  append_val(v_num_reads,             c_num_reads);
  append_val(v_uncompressed_names_sz, c_uncompressed_names_sz);
  append_val(v_uncompressed_seqs_sz,  c_uncompressed_seqs_sz);
  append_val(v_name_prefix_sz,        c_name_prefix_sz);
  append_val(v_cname_lengths_sz,      c_cname_lengths_sz);
  append_val(v_cnames_sz,             c_cnames_sz);
  append_val(v_cseq_lengths_sz,       c_cseq_lengths_sz);
  append_val(v_cseqs_sz,              c_cseqs_sz);
  append_val(v_cquals_sz,             c_cquals_sz);
  
  return out;
}


std::array<uint8_t, 16> generate_file_header() {
  constexpr uint8_t FRAQ_FORMAT_VERSION = 1;
  bool be = is_big_endian();
  // Write a file header with magic bits and version
  std::array<uint8_t, 16> file_header = 
    {'F','R','A','Q', // magic bits
     0,0,0,0,
     0,0,0,0,
     0,0,
     // big-endian flag, files written in big-endian cant be read on little-endian systems and vice versa
     static_cast<uint8_t>(be),
     FRAQ_FORMAT_VERSION}; // version 1
  return file_header;
}

void write_fraq(std::vector<Read> &&reads, std::ofstream &ofs, const int compress_level) {
  auto file_header = generate_file_header(); // std::array<char, ...>
  ofs.write(reinterpret_cast<char*>(file_header.data()), file_header.size());
  
  size_t index = 0;
  while (index < reads.size()) {
    size_t n = std::min<size_t>(reads.size() - index, MAX_BLOCKSIZE);
    CompressedReadBlock block = compress_reads(reads.data() + index, n, compress_level);
    
    std::vector<char> header = generate_block_header(block);
    ofs.write(header.data(), header.size());
    
    ofs.write(block.name_prefix.data(), block.name_prefix.size());
    ofs.write(block.compressed_name_lengths.data(), block.compressed_name_lengths.size());
    ofs.write(block.compressed_names.data(), block.compressed_names.size());
    ofs.write(block.compressed_seq_lengths.data(), block.compressed_seq_lengths.size());
    ofs.write(block.compressed_seqs.data(), block.compressed_seqs.size());
    ofs.write(block.compressed_quals.data(), block.compressed_quals.size());
    
    index += n;
  }
}

////////////////////////////////////////////////////////////////////////////////
// read

inline uint8_t read_fraq_header(std::istream& is, uint8_t supported_version) {
  std::array<uint8_t, 16> h{};
  if (!is.read(reinterpret_cast<char*>(h.data()), h.size())) {
    throw std::runtime_error("failed to read FRAQ header");
  }
  
  if (!(h[0]=='F' && h[1]=='R' && h[2]=='A' && h[3]=='Q')) {
    throw std::runtime_error("bad FRAQ magic");
  }
  
  if (h[14] != static_cast<uint8_t>(is_big_endian())) {
    throw std::runtime_error("endianness mismatch (v1 not portable)");
  }
  
  if (h[15] == 0 || h[15] > supported_version) {
    throw std::runtime_error("unsupported FRAQ version");
  }
  
  return h[15];
}


CompressedReadBlock read_block(std::istream& is) {
  using ES = DynArray::ElemSize;
  
  auto byte_count = [](ES es) -> size_t {
    return es == ES::B8 ? 1 : es == ES::B16 ? 2 : es == ES::B32 ? 4 : 8;
  };
  
  auto read_uint = [&](ES es) -> uint64_t {
    uint64_t v = 0;
    const size_t n = byte_count(es);
    if (!is.read(reinterpret_cast<char*>(&v), n)) {
      throw std::runtime_error("short read while reading header value");
    }
    return v; // native-endian
  };
  
  auto rd_exact = [&](void* p, size_t n) {
    if (!is.read(static_cast<char*>(p), n)) {
      throw std::runtime_error("short read while reading block payload");
    }
  };
  
  // meta
  uint32_t meta = 0;
  if (!is.read(reinterpret_cast<char*>(&meta), 4)) {
    throw std::runtime_error("short read while reading block meta");
  }
  
  auto cls = [&](unsigned shift) -> ES {
    return static_cast<ES>((meta >> shift) & 0x3u);
  };
  
  const ES c_num_reads             = cls(0);
  const ES c_uncompressed_names_sz = cls(2);
  const ES c_uncompressed_seqs_sz  = cls(4);
  const ES c_name_prefix_sz        = cls(6);
  const ES c_cname_lengths_sz      = cls(8);
  const ES c_cnames_sz             = cls(10);
  const ES c_cseq_lengths_sz       = cls(12);
  const ES c_cseqs_sz              = cls(14);
  const ES c_cquals_sz             = cls(16);
  
  const ES name_len_elem_size = cls(18);
  const ES seq_len_elem_size  = cls(20);
  const bool use_bit_pack     = ((meta >> 22) & 0x1u) != 0;
  
  const uint64_t v_num_reads             = read_uint(c_num_reads);
  const uint64_t v_uncompressed_names_sz = read_uint(c_uncompressed_names_sz);
  const uint64_t v_uncompressed_seqs_sz  = read_uint(c_uncompressed_seqs_sz);
  const uint64_t v_name_prefix_sz        = read_uint(c_name_prefix_sz);
  const uint64_t v_cname_lengths_sz      = read_uint(c_cname_lengths_sz);
  const uint64_t v_cnames_sz             = read_uint(c_cnames_sz);
  const uint64_t v_cseq_lengths_sz       = read_uint(c_cseq_lengths_sz);
  const uint64_t v_cseqs_sz              = read_uint(c_cseqs_sz);
  const uint64_t v_cquals_sz             = read_uint(c_cquals_sz);
  
  CompressedReadBlock b;
  b.name_len_elem_size      = name_len_elem_size;
  b.seq_len_elem_size       = seq_len_elem_size;
  b.use_bit_pack            = use_bit_pack;
  b.num_reads               = static_cast<uint16_t>(v_num_reads); // writer guarantees <= 65535
  b.uncompressed_names_size = v_uncompressed_names_sz;
  b.uncompressed_seqs_size  = v_uncompressed_seqs_sz;
  
  b.name_prefix.resize(static_cast<size_t>(v_name_prefix_sz));
  rd_exact(b.name_prefix.data(), b.name_prefix.size());
  
  b.compressed_name_lengths.resize(static_cast<size_t>(v_cname_lengths_sz));
  rd_exact(b.compressed_name_lengths.data(), b.compressed_name_lengths.size());
  
  b.compressed_names.resize(static_cast<size_t>(v_cnames_sz));
  rd_exact(b.compressed_names.data(), b.compressed_names.size());
  
  b.compressed_seq_lengths.resize(static_cast<size_t>(v_cseq_lengths_sz));
  rd_exact(b.compressed_seq_lengths.data(), b.compressed_seq_lengths.size());
  
  b.compressed_seqs.resize(static_cast<size_t>(v_cseqs_sz));
  rd_exact(b.compressed_seqs.data(), b.compressed_seqs.size());
  
  b.compressed_quals.resize(static_cast<size_t>(v_cquals_sz));
  rd_exact(b.compressed_quals.data(), b.compressed_quals.size());
  
  return b;
}

std::vector<Read> decompress_block(const CompressedReadBlock& b) {
  using ES = DynArray::ElemSize;
  
  auto bytes_per = [](ES es) -> size_t {
    return es == ES::B8 ? 1 : es == ES::B16 ? 2 : es == ES::B32 ? 4 : 8;
  };
  auto load_u64 = [&](const std::vector<char>& buf, size_t idx, ES es) -> uint64_t {
    uint64_t v = 0;
    const size_t stride = bytes_per(es);
    std::memcpy(&v, buf.data() + idx * stride, stride); // native-endian
    return v;
  };
  auto zstd_decompress_exact = [](std::vector<char>& out, size_t out_size,
                                  const std::vector<char>& in) {
    out.resize(out_size);
    size_t got = ZSTD_decompress(out.data(), out.size(), in.data(), in.size());
    if (ZSTD_isError(got) || got != out_size) {
      throw std::runtime_error("zstd decompress failed or size mismatch");
    }
  };
  
  const size_t n_reads = b.num_reads;
  
  // 1) Decompress length arrays first (we need them to size quals and rebuild fields)
  const size_t name_len_bytes = bytes_per(b.name_len_elem_size) * n_reads;
  const size_t seq_len_bytes  = bytes_per(b.seq_len_elem_size)  * n_reads;
  
  std::vector<char> name_len_buf;
  std::vector<char> seq_len_buf;
  zstd_decompress_exact(name_len_buf, name_len_bytes, b.compressed_name_lengths);
  zstd_decompress_exact(seq_len_buf,  seq_len_bytes,  b.compressed_seq_lengths);
  
  // Compute totals for validation and for quals sizing
  uint64_t total_name_tail = 0;
  uint64_t total_seq_len   = 0;
  for (size_t i = 0; i < n_reads; ++i) {
    total_name_tail += load_u64(name_len_buf, i, b.name_len_elem_size);
    total_seq_len   += load_u64(seq_len_buf,  i, b.seq_len_elem_size);
  }
  
  // Optional sanity checks (kept simple)
  if (total_name_tail != b.uncompressed_names_size) {
    throw std::runtime_error("names payload size mismatch");
  }
  
  // 2) Decompress payloads
  std::vector<char> names_payload;
  std::vector<char> seqs_payload;
  std::vector<char> quals_payload;
  
  zstd_decompress_exact(names_payload, static_cast<size_t>(b.uncompressed_names_size),
                        b.compressed_names);
  zstd_decompress_exact(seqs_payload,  static_cast<size_t>(b.uncompressed_seqs_size),
                        b.compressed_seqs);
  zstd_decompress_exact(quals_payload, static_cast<size_t>(total_seq_len),
                        b.compressed_quals);
  
  // 3) Rebuild reads
  std::vector<Read> out;
  out.reserve(n_reads);
  
  const char* pname = names_payload.data();
  const char* pqual = quals_payload.data();
  
  // inverse code map for bit-packed sequences
  static const char inv_code[16] = {
    'A','C','G','T','R','Y','S','W','K','M','B','D','H','V','N','U'
  };
  const unsigned char* pseq_packed = reinterpret_cast<const unsigned char*>(seqs_payload.data());
  const char* pseq_raw = seqs_payload.data();
  uint64_t pos_nibbles = 0; // only used if bit-packed
  
  for (size_t i = 0; i < n_reads; ++i) {
    const uint64_t name_tail_len = load_u64(name_len_buf, i, b.name_len_elem_size);
    const uint64_t seq_len       = load_u64(seq_len_buf,  i, b.seq_len_elem_size);
    
    Read r;
    r.name.reserve(b.name_prefix.size() + static_cast<size_t>(name_tail_len));
    r.name = b.name_prefix;
    r.name.append(pname, static_cast<size_t>(name_tail_len));
    pname += name_tail_len;
    
    r.seq.resize(static_cast<size_t>(seq_len));
    if (b.use_bit_pack) {
      for (size_t k = 0; k < seq_len; ++k, ++pos_nibbles) {
        const size_t j = static_cast<size_t>(pos_nibbles >> 1);
        const uint8_t byte = pseq_packed[j];
        const uint8_t code = (pos_nibbles & 1) ? (byte & 0x0Fu) : ((byte >> 4) & 0x0Fu);
        r.seq[k] = inv_code[code];
      }
    } else {
      r.seq.assign(pseq_raw, pseq_raw + seq_len);
      pseq_raw += seq_len;
    }
    
    r.qual.assign(pqual, pqual + seq_len);
    pqual += seq_len;
    
    out.push_back(std::move(r));
  }
  return out;
}

std::vector<Read> read_fraq(std::istream& is) {
  // validate file header (native-endian v1)
  constexpr uint8_t SUPPORTED_VERSION = 1;
  (void)read_fraq_header(is, SUPPORTED_VERSION);
  
  std::vector<Read> all;
  
  while(true) {
    // stop cleanly at EOF between blocks
    int next = is.peek();
    if (next == std::char_traits<char>::eof()) break;
    
    CompressedReadBlock blk = read_block(is);
    std::vector<Read> part = decompress_block(blk);
    
    // append (move) reads from this block
    all.insert(all.end(),
               std::make_move_iterator(part.begin()),
               std::make_move_iterator(part.end()));
  }
  
  return all;
}

} // namespace fraq_internal

#endif // FRAQ_FORMAT_H