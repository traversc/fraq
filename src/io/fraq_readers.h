#ifndef __FRAQ_READERS_H
#define __FRAQ_READERS_H

#include <cstddef>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cstdio>
#include <array>
#include <zlib.h>
#include <zstd.h>
#include <unordered_map>
#include "fraq_defines.h"
#include "fraq_format.h"

#if FRAQ_HAS_FIFO
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#endif

namespace fraq_internal {

using fraq::ends_with;

static constexpr size_t DEFAULT_BUF_SIZE = 64 * 1024;
static constexpr std::array<uint8_t, 2> GZIP_MAGIC{{0x1f, 0x8b}};
static constexpr std::array<uint8_t, 4> ZSTD_MAGIC{{0x28, 0xb5, 0x2f, 0xfd}};

template <size_t N>
inline void ensure_magic_prefix(const std::string& path,
                                const std::array<uint8_t, N>& magic,
                                const char* format_name) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error(std::string("Failed to open ") + format_name + " file: " + path);
    }
    std::array<char, N> buf{};
    in.read(buf.data(), static_cast<std::streamsize>(N));
    if (in.gcount() < static_cast<std::streamsize>(N)) {
        throw std::runtime_error(std::string("File too short to be valid ") + format_name + ": " + path);
    }
    for (size_t i = 0; i < N; ++i) {
        if (static_cast<uint8_t>(buf[i]) != magic[i]) {
            throw std::runtime_error(std::string("File does not appear to be valid ") + format_name + ": " + path);
        }
    }
}

struct IReader {
    virtual ~IReader() = default;

    // Byte/line primitives for FASTQ-based readers.
    virtual size_t read(char* dst, size_t cnt) = 0;
    virtual bool   readLine(std::string& out) = 0;

    // Unified API: produce up to block_size reads.
    // Default implementation parses FASTQ lines via readLine().
    virtual void read_block(std::vector<fraq::Read>& block, size_t block_size) {
        block.clear();
        block.reserve(block_size);
        while (block.size() < block_size) {
            fraq::Read r;
            if (!readLine(r.name)) break;
            if (r.name.empty() || r.name[0] != '@') {
                throw std::runtime_error("Invalid FASTQ read name: " + r.name);
            }
            r.name.erase(0, 1);
            if (!readLine(r.seq)) break;
            std::string plus_line;
            if (!readLine(plus_line)) break;
            if (!readLine(r.qual)) break;
            block.push_back(std::move(r));
        }
    }
};

// gzip reader
struct GzipReader : IReader {
    gzFile             _gz;
    std::vector<char>  _buf;
    size_t             _pos;
    size_t             _len;

    explicit GzipReader(const std::string& path, size_t bufSize)
        : _gz(nullptr),
          _buf(bufSize),
          _pos(0),
          _len(0)
    {
        ensure_magic_prefix(path, GZIP_MAGIC, "gzip");
        _gz = gzopen(path.c_str(), "rb");
        if (!_gz) throw std::runtime_error("gzip open failed: " + path);
    }

    ~GzipReader() override {
        if (_gz) gzclose(_gz);
    }

    size_t read(char* dst, size_t cnt) override {
        size_t tot = 0;
        while (tot < cnt) {
            if (_pos == _len) {
                int n = gzread(_gz, _buf.data(),
                               static_cast<unsigned>(_buf.size()));
                if (n <= 0) break;
                _pos = 0; _len = static_cast<size_t>(n);
            }
            size_t chunk = std::min(cnt - tot, _len - _pos);
            std::memcpy(dst + tot, _buf.data() + _pos, chunk);
            tot  += chunk;
            _pos += chunk;
        }
        return tot;
    }

    bool readLine(std::string& out) override {
        out.clear();
        while (true) {
            if (_pos == _len) {
                int n = gzread(_gz, _buf.data(),
                               static_cast<unsigned>(_buf.size()));
                if (n <= 0) return !out.empty();
                _pos = 0; _len = static_cast<size_t>(n);
            }
            char*  start = _buf.data() + _pos;
            size_t rem   = _len - _pos;
            void*  nl    = std::memchr(start, '\n', rem);
            if (nl) {
                size_t chunk = static_cast<char*>(nl) - start;
                out.append(start, chunk);
                _pos += chunk + 1;
                return true;
            } else {
                out.append(start, rem);
                _pos += rem;
            }
        }
    }
};

// zstd reader
struct ZstdReader : IReader {
    FILE*              _f;
    ZSTD_DCtx*         _dctx;
    std::vector<char>  _inBuf;
    std::vector<char>  _outBuf;
    size_t             _inPos;
    size_t             _inSize;
    size_t             _outPos;
    size_t             _outSize;

    explicit ZstdReader(const std::string& path, size_t bufSize)
        : _f(nullptr),
          _dctx(nullptr),
          _inBuf(bufSize),
          _outBuf(bufSize),
          _inPos(0),
          _inSize(0),
          _outPos(0),
          _outSize(0)
    {
        ensure_magic_prefix(path, ZSTD_MAGIC, "zstd");
        _f = std::fopen(path.c_str(), "rb");
        if (!_f) throw std::runtime_error("zstd open failed: " + path);

        _dctx = ZSTD_createDCtx();
        if (!_dctx) {
            std::fclose(_f);
            throw std::runtime_error("ZSTD_createDCtx failed");
        }
        size_t init = ZSTD_initDStream(_dctx);
        if (ZSTD_isError(init)) {
            ZSTD_freeDCtx(_dctx);
            std::fclose(_f);
            throw std::runtime_error("ZSTD_initDStream failed: " +
                                     std::string(ZSTD_getErrorName(init)));
        }
    }

    ~ZstdReader() override {
        if (_dctx) ZSTD_freeDCtx(_dctx);
        if (_f)    std::fclose(_f);
    }

    size_t read(char* dst, size_t cnt) override {
        size_t tot = 0;
        while (tot < cnt) {
            if (_outPos == _outSize) {
                if (_inPos == _inSize) {
                    _inSize = std::fread(_inBuf.data(), 1,
                                         _inBuf.size(), _f);
                    _inPos = 0;
                    if (_inSize == 0) break;
                }
                ZSTD_inBuffer  in  = {_inBuf.data() + _inPos, _inSize - _inPos, 0};
                ZSTD_outBuffer out = {_outBuf.data(), _outBuf.size(), 0};
                size_t ret = ZSTD_decompressStream(_dctx, &out, &in);
                if (ZSTD_isError(ret))
                    throw std::runtime_error("ZSTD_decompressStream failed: " +
                                             std::string(ZSTD_getErrorName(ret)));
                _inPos   += in.pos;
                _outSize = out.pos;
                _outPos   = 0;
                if (_outSize == 0) continue;
            }
            size_t chunk = std::min(cnt - tot, _outSize - _outPos);
            std::memcpy(dst + tot, _outBuf.data() + _outPos, chunk);
            tot     += chunk;
            _outPos += chunk;
        }
        return tot;
    }

    bool readLine(std::string& out) override {
        out.clear();
        char c;
        while (read(&c, 1) == 1) {
            if (c == '\n') return true;
            out.push_back(c);
        }
        return !out.empty();
    }
};

// plain-text reader
struct PlainReader : IReader {
    std::ifstream _ifs;

    explicit PlainReader(const std::string& path)
      : _ifs(path, std::ios::binary)
    {
        if (!_ifs) throw std::runtime_error("plain open failed: " + path);
    }

    ~PlainReader() override = default;

    size_t read(char* dst, size_t cnt) override {
        _ifs.read(dst, cnt);
        return static_cast<size_t>(_ifs.gcount());
    }

    bool readLine(std::string& out) override {
        return static_cast<bool>(std::getline(_ifs, out));
    }
};

#if FRAQ_HAS_FIFO
// FIFO reader for plain FASTQ
struct FifoReader : IReader {
    int                _fd;
    std::vector<char>  _buf;
    size_t             _pos;
    size_t             _len;

    explicit FifoReader(const std::string& path, size_t bufSize)
        : _fd(-1),
          _buf(bufSize),
          _pos(0),
          _len(0)
    {
        _fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (_fd < 0) throw std::runtime_error("fifo open failed: " + path);
    }

    ~FifoReader() override {
        if (_fd >= 0) ::close(_fd);
    }

    size_t read(char* dst, size_t cnt) override {
        size_t tot = 0;
        while (tot < cnt) {
            if (_pos == _len) {
                ssize_t n;
                do {
                    n = ::read(_fd, _buf.data(), static_cast<unsigned>(_buf.size()));
                } while (n < 0 && errno == EINTR);
                if (n <= 0) break;
                _pos = 0; _len = static_cast<size_t>(n);
            }
            const size_t chunk = std::min(cnt - tot, _len - _pos);
            std::memcpy(dst + tot, _buf.data() + _pos, chunk);
            tot  += chunk;
            _pos += chunk;
        }
        return tot;
    }

    bool readLine(std::string& out) override {
        out.clear();
        while (true) {
            if (_pos == _len) {
                ssize_t n;
                do {
                    n = ::read(_fd, _buf.data(), static_cast<unsigned>(_buf.size()));
                } while (n < 0 && errno == EINTR);
                if (n <= 0) return !out.empty();
                _pos = 0; _len = static_cast<size_t>(n);
            }
            char*  start = _buf.data() + _pos;
            size_t rem   = _len - _pos;
            void*  nl    = std::memchr(start, '\n', rem);
            if (nl) {
                size_t chunk = static_cast<char*>(nl) - start;
                out.append(start, chunk);
                _pos += chunk + 1;
                return true;
            } else {
                out.append(start, rem);
                _pos += rem;
            }
        }
    }
};
#else
struct FifoReader : IReader {
    explicit FifoReader(const std::string& path, size_t) {
        (void)path;
        throw std::runtime_error("FIFO input is not supported on Windows");
    }
    size_t read(char*, size_t) override { return 0; }
    bool readLine(std::string&) override { return false; }
};
#endif

// FRAQ file reader (block-oriented)
struct FraqfFileReader : IReader {
    std::ifstream               _ifs;
    std::vector<fraq::Read>     _buf;
    size_t                      _pos;

    explicit FraqfFileReader(const std::string& path)
      : _ifs(path, std::ios::binary),
        _buf(),
        _pos(0)
    {
        if (!_ifs) throw std::runtime_error("fraq open failed: " + path);
        static const uint8_t SUPPORTED_VERSION = 1;
        (void)fraq_internal::read_fraq_header(_ifs, SUPPORTED_VERSION);
    }

    size_t read(char* /*dst*/, size_t /*cnt*/) override {
        throw std::runtime_error("read() unsupported for FRAQ; use read_block()");
    }
    bool readLine(std::string& /*out*/) override {
        return false;
    }

    void read_block(std::vector<fraq::Read>& out, size_t block_size) override {
        out.clear();
        out.reserve(block_size);

        while (out.size() < block_size) {
            if (_pos == _buf.size()) {
                int next = _ifs.peek();
                if (next == std::char_traits<char>::eof()) break;

                CompressedReadBlock blk = fraq_internal::read_block(_ifs);
                _buf = fraq_internal::decompress_block(blk);
                _pos = 0;
            }

            const size_t can = std::min(block_size - out.size(),
                                        _buf.size() - _pos);

            out.insert(out.end(),
                       std::make_move_iterator(_buf.begin() + _pos),
                       std::make_move_iterator(_buf.begin() + _pos + can));

            _pos += can;
        }
    }
};

// In-memory FRAQ reader backed by a concurrent map.
// Assumes the mapped vector for the mem_key remains present and is not modified
// in a way that replaces the vector while this reader is active.
struct FraqfMemReader : IReader {
    const std::vector<CompressedReadBlock>* _blocks;
    size_t                                  _blk_idx;
    std::vector<fraq::Read>                 _spill;
    size_t                                  _spill_pos;

    FraqfMemReader(
        std::unordered_map<std::string, std::vector<CompressedReadBlock>>& store,
        const std::string& mem_key
    )
        : _blocks(nullptr),
          _blk_idx(0),
          _spill(),
          _spill_pos(0)
    {
        if (!ends_with(mem_key, ".mem")) {
            throw std::runtime_error("FraqfMemReader key must end with .mem: " + mem_key);
        }
        auto it = store.find(mem_key);
        if (it == store.end()) {
            throw std::runtime_error("FraqfMemReader key not found: " + mem_key);
        }
        // Cache a pointer to the mapped vector (see discussion about stability).
        _blocks = &it->second;
    }

    size_t read(char* /*dst*/, size_t /*cnt*/) override {
        throw std::runtime_error("read() unsupported for FraqfMemReader; use read_block()");
    }
    bool readLine(std::string& /*out*/) override {
        return false;
    }

    void read_block(std::vector<fraq::Read>& out, size_t block_size) override {
        out.clear();
        out.reserve(block_size);

        // Drain any spill first
        while (_spill_pos < _spill.size() && out.size() < block_size) {
            out.push_back(std::move(_spill[_spill_pos++]));
        }
        if (_spill_pos == _spill.size()) {
            _spill.clear();
            _spill_pos = 0;
        }

        while (out.size() < block_size) {
            if (_blk_idx >= _blocks->size()) break;

            const CompressedReadBlock& blk = (*_blocks)[_blk_idx++];
            std::vector<fraq::Read> decompressed = fraq_internal::decompress_block(blk);

            const size_t need = block_size - out.size();
            const size_t take = std::min(need, decompressed.size());

            // Move the first 'take' reads to out
            out.insert(out.end(),
                       std::make_move_iterator(decompressed.begin()),
                       std::make_move_iterator(decompressed.begin() + take));

            // Keep any leftovers as spill
            if (take < decompressed.size()) {
                _spill = std::move(decompressed);
                _spill_pos = take;
                break;
            }
        }
    }
};

// runtime dispatcher
struct MultiReader {
    std::unique_ptr<IReader> impl;

    explicit MultiReader(const std::string& path) {
        if      (ends_with(path, ".fraq")) impl = std::make_unique<FraqfFileReader>(path);
        else if (ends_with(path, ".gz"))   impl = std::make_unique<GzipReader>(path, DEFAULT_BUF_SIZE);
        else if (ends_with(path, ".zst"))  impl = std::make_unique<ZstdReader>(path, DEFAULT_BUF_SIZE);
#if FRAQ_HAS_FIFO
        else if (ends_with(path, ".fifo")) impl = std::make_unique<FifoReader>(path, DEFAULT_BUF_SIZE);
#else
        else if (ends_with(path, ".fifo")) {
            throw std::runtime_error("FIFO inputs are unavailable on this platform");
        }
#endif
        else                               impl = std::make_unique<PlainReader>(path);
    }

    // Overload for in-memory FRAQ store
    MultiReader(
        const std::string& mem_key,
        std::unordered_map<std::string, std::vector<CompressedReadBlock>>& store
    ) {
        impl = std::make_unique<FraqfMemReader>(store, mem_key);
    }

    void read_block(std::vector<fraq::Read>& block, size_t block_size) {
        impl->read_block(block, block_size);
    }
};

} // namespace fraq_internal

#endif // __FRAQ_READERS_H
