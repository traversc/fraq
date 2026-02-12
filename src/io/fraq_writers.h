#ifndef __FRAQ_WRITERS_H
#define __FRAQ_WRITERS_H

#include <cstddef>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cstdio>
#include <cerrno>
#include <zlib.h>
#include <zstd.h>
#include <unordered_map>
#include <deque>
#include "fraq_defines.h"
#include "fraq_format.h"

// FIFO compile time support
#if FRAQ_HAS_FIFO
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <poll.h>
#endif

namespace fraq_internal {
using fraq::ends_with;

// writer interface
struct IWriter {
    virtual ~IWriter() = default;
    virtual void write(const char *data, size_t len) = 0;
    virtual void writeLine(const std::string &line) = 0;
    virtual void flush() = 0;
    virtual bool has_pending() const { return false; }
    virtual void close() {}
};

// gzip writer
struct GzipWriter : IWriter {
    gzFile            _gz;
    std::vector<char> _buf;
    size_t            _bufSize;

    explicit GzipWriter(const std::string &path,
                        int compress_level,
                        size_t bufSize = 64*1024)
      : _gz(gzopen(path.c_str(), ("wb" + std::to_string(compress_level)).c_str()))
      , _buf()
      , _bufSize(bufSize)
    {
        if (!_gz) throw std::runtime_error("Failed to open " + path);
        _buf.reserve(_bufSize);
    }

    ~GzipWriter() override {
        try { close(); } catch(...) {}
    }

    void write(const char *data, size_t len) override {
        size_t done = 0;
        while (done < len) {
            size_t space = _bufSize - _buf.size();
            size_t chunk = std::min(space, len - done);
            _buf.insert(_buf.end(), data + done, data + done + chunk);
            done += chunk;
            if (_buf.size() == _bufSize) flush();
        }
    }

    void writeLine(const std::string &line) override {
        write(line.data(), line.size());
        write("\n", 1);
    }

    void flush() override {
        if (_buf.empty()) return;
        if (gzwrite(_gz, _buf.data(), _buf.size()) <= 0)
            throw std::runtime_error("gzwrite failed");
        _buf.clear();
        gzflush(_gz, Z_SYNC_FLUSH);
    }

    void close() override {
        if (!_gz) return;
        flush();
        gzclose(_gz);
        _gz = nullptr;
    }
};

// zstd writer
struct ZstdWriter : IWriter {
    FILE*              _f;
    ZSTD_CCtx*         _cctx;
    std::vector<char>  _inBuf;
    std::vector<char>  _outBuf;
    size_t             _bufSize;
    size_t             _inPos{0};

    explicit ZstdWriter(const std::string &path,
                        int compress_level,
                        size_t bufSize = 64*1024)
      : _f(std::fopen(path.c_str(), "wb"))
      , _cctx(ZSTD_createCCtx())
      , _inBuf(bufSize)
      , _outBuf(ZSTD_CStreamOutSize())
      , _bufSize(bufSize)
    {
        if (!_f) {
            throw std::runtime_error("Failed to open " + path);
        }
        if (!_cctx) {
            std::fclose(_f);
            throw std::runtime_error("ZSTD_createCCtx failed");
        }
        size_t init = ZSTD_initCStream(_cctx, compress_level);
        if (ZSTD_isError(init)) {
            ZSTD_freeCCtx(_cctx);
            std::fclose(_f);
            throw std::runtime_error(std::string("ZSTD_initCStream failed: ")
                                     + ZSTD_getErrorName(init));
        }
    }

    ~ZstdWriter() override {
        try { close(); } catch(...) {}
    }

    void write(const char* data, size_t len) override {
        size_t done = 0;
        while (done < len) {
            size_t space = _bufSize - _inPos;
            size_t chunk = std::min(space, len - done);
            std::memcpy(_inBuf.data() + _inPos, data + done, chunk);
            _inPos += chunk;
            done   += chunk;
            if (_inPos == _bufSize) {
                compressChunk(); // compress full input buffer
            }
        }
    }

    void writeLine(const std::string &line) override {
        write(line.data(), line.size());
        static const char nl = '\n';
        write(&nl, 1);
    }

    void flush() override {
        if (_inPos > 0) {
            compressChunk(); // push remaining input
        }

        // finish the frame
        ZSTD_outBuffer out = { _outBuf.data(), _outBuf.size(), 0 };
        size_t ret;
        do {
            ret = ZSTD_endStream(_cctx, &out);
            if (ZSTD_isError(ret)) {
                throw std::runtime_error(std::string("ZSTD_endStream failed: ")
                                         + ZSTD_getErrorName(ret));
            }
            if (out.pos > 0) {
                size_t written = std::fwrite(_outBuf.data(), 1, out.pos, _f);
                if (written != out.pos) {
                    throw std::runtime_error("fwrite failed");
                }
                out.pos = 0;
            }
        } while (ret > 0);

        std::fflush(_f);
    }

    void close() override {
        if (!_cctx && !_f) return;
        flush();
        if (_cctx) {
            ZSTD_freeCCtx(_cctx);
            _cctx = nullptr;
        }
        if (_f) {
            std::fclose(_f);
            _f = nullptr;
        }
    }

private:
    void compressChunk() {
        ZSTD_inBuffer  in  = { _inBuf.data(), _inPos, 0 };
        ZSTD_outBuffer out = { _outBuf.data(), _outBuf.size(), 0 };

        while (in.pos < in.size) {
            size_t ret = ZSTD_compressStream(_cctx, &out, &in);
            if (ZSTD_isError(ret)) {
                throw std::runtime_error(std::string("ZSTD_compressStream failed: ")
                                         + ZSTD_getErrorName(ret));
            }
            if (out.pos > 0) {
                size_t written = std::fwrite(_outBuf.data(), 1, out.pos, _f);
                if (written != out.pos) {
                    throw std::runtime_error("fwrite failed");
                }
                out.pos = 0; // reset for next iteration
            }
        }
        _inPos = 0; // input buffer fully consumed
    }
};

// plain-text writer
struct PlainWriter : IWriter {
    std::ofstream _ofs;

    explicit PlainWriter(const std::string &path)
      : _ofs(path, std::ios::binary)
    {
        if (!_ofs) throw std::runtime_error("Failed to open " + path);
    }

    ~PlainWriter() override {
        try { close(); } catch(...) {}
    }

    void write(const char *data, size_t len) override {
        _ofs.write(data, len);
        if (!_ofs) throw std::runtime_error("write failed");
    }

    void writeLine(const std::string &line) override {
        write(line.data(), line.size());
        write("\n", 1);
    }

    void flush() override {
        _ofs.flush();
        if (!_ofs) throw std::runtime_error("flush failed");
    }

    void close() override {
        flush();
        _ofs.close();
    }
};

#if FRAQ_HAS_FIFO
struct FifoWriter : IWriter {
    std::string path;
    int fd{-1};
    std::deque<std::vector<char>> pending;
    size_t front_offset{0};

    explicit FifoWriter(const std::string& fifo_path) : path(fifo_path) {
        struct stat st{};
        if (::lstat(path.c_str(), &st) == 0) {
            if (!S_ISFIFO(st.st_mode)) {
                throw std::runtime_error("Path exists and is not a FIFO: " + path);
            }
        } else {
            if (::mkfifo(path.c_str(), 0666) != 0) {
                throw std::runtime_error("mkfifo failed for " + path + ": errno " + std::to_string(errno));
            }
        }

        // Best-effort open; if no reader is attached yet we will retry later.
        try_open_nonblocking();
    }

    void write(const char* data, size_t len) override {
        if (len == 0) return;
        if (!pending.empty()) {
            buffer_append(data, len);
            flush();
            return;
        }
        size_t off = 0;
        while (off < len) {
            if (fd == -1 && !try_open_nonblocking()) {
                buffer_append(data + off, len - off);
                return;
            }
            ssize_t n = ::write(fd, data + off, len - off);
            if (n > 0) {
                off += static_cast<size_t>(n);
                continue;
            }
            if (n == -1 && errno == EINTR) continue;
            if (n == -1 && (errno == EAGAIN || errno == EWOULDBLOCK)) {
                buffer_append(data + off, len - off);
                return;
            }
            if (n == -1 && errno == EPIPE) {
                close_fd();
                throw std::runtime_error("FIFO reader closed unexpectedly while writing");
            }
            throw std::runtime_error("write to FIFO failed: errno " + std::to_string(errno));
        }
    }

    void writeLine(const std::string& line) override {
        write(line.data(), line.size());
        static const char nl = '\n';
        write(&nl, 1);
    }

    bool has_pending() const override {
        return !pending.empty();
    }

    void flush() override {
        if (pending.empty()) {
            // Still attempt to establish a reader for future writes.
            try_open_nonblocking();
            return;
        }
        if (!try_open_nonblocking()) {
            return; // nothing to do until a reader attaches
        }
        while (!pending.empty()) {
            std::vector<char>& chunk = pending.front();
            const size_t remaining = chunk.size() - front_offset;
            ssize_t n = ::write(fd, chunk.data() + front_offset, remaining);
            if (n > 0) {
                front_offset += static_cast<size_t>(n);
                if (front_offset == chunk.size()) {
                    pending.pop_front();
                    front_offset = 0;
                }
                continue;
            }
            if (n == -1 && errno == EINTR) continue;
            if (n == -1 && (errno == EAGAIN || errno == EWOULDBLOCK)) {
                return;
            }
            if (n == -1 && errno == EPIPE) {
                close_fd();
                throw std::runtime_error("FIFO reader closed unexpectedly while flushing");
            }
            throw std::runtime_error("write to FIFO failed: errno " + std::to_string(errno));
        }
        front_offset = 0;
    }

    void close() override {
        if (!pending.empty()) {
            throw std::runtime_error("Attempted to close FIFO writer with pending data");
        }
        close_fd();
    }
  private:
    bool try_open_nonblocking() {
        if (fd != -1) return true;
        int open_flags = O_WRONLY | O_NONBLOCK;
#ifdef O_CLOEXEC
        open_flags |= O_CLOEXEC;
#endif
        while (true) {
            int new_fd = ::open(path.c_str(), open_flags);
            if (new_fd != -1) {
                fd = new_fd;
                return true;
            }
            if (errno == EINTR) continue;
            if (errno == ENXIO) {
                return false; // reader not ready yet
            }
            throw std::runtime_error("open FIFO for write failed: errno " + std::to_string(errno));
        }
    }

    void close_fd() noexcept {
        if (fd != -1) {
            ::close(fd);
            fd = -1;
        }
    }

    void buffer_append(const char* data, size_t len) {
        if (len == 0) return;
        pending.emplace_back(data, data + len);
        if (pending.size() == 1) {
            front_offset = 0;
        }
    }
};
#else
struct FifoWriter : IWriter {
    explicit FifoWriter(const std::string& path) {
        (void)path;
        throw std::runtime_error("FIFO output is not supported on Windows");
    }
    void write(const char*, size_t) override {}
    void writeLine(const std::string&) override {}
    bool has_pending() const override { return false; }
    void flush() override {}
    void close() override {}
};
#endif


// dispatcher that picks writer by extension
struct MultiWriter {
    std::unique_ptr<IWriter> impl;

    explicit MultiWriter(const std::string &path, int compress_level) {
        if      (ends_with(path, ".gz"))   impl = std::make_unique<GzipWriter>(path, compress_level);
        else if (ends_with(path, ".zst"))  impl = std::make_unique<ZstdWriter>(path, compress_level);
#if FRAQ_HAS_FIFO
        else if (ends_with(path, ".fifo")) impl = std::make_unique<FifoWriter>(path);
#else
        else if (ends_with(path, ".fifo")) {
            throw std::runtime_error("FIFO outputs are unavailable on this platform");
        }
#endif
        else                               impl = std::make_unique<PlainWriter>(path);
    }
    bool has_pending() const {
        return impl && impl->has_pending();
    }
    void flush() {
        impl->flush();
    }
    void close() {
        if (impl) impl->close();
    }
    // add this method inside struct MultiWriter
    void write_block(const std::vector<fraq::Read>& reads) {
        for (const auto& r : reads) {
            impl->write("@", 1);
            impl->writeLine(r.name);
            impl->writeLine(r.seq);
            impl->writeLine("+");
            impl->writeLine(r.qual);
        }
    }
};


// FRAQ writer interface
struct IFraqfWriter {
    virtual ~IFraqfWriter() = default;
    virtual void write_compressed_block(CompressedReadBlock&& block) = 0; // destructive, must std::move
    virtual void flush() = 0;
    virtual void close() = 0;
};

// .fraq on-disk writer
struct FraqfFileWriter : IFraqfWriter {
    std::ofstream ofs;

    explicit FraqfFileWriter(const std::string& path)
        : ofs(path, std::ios::binary)
    {
        if (!ends_with(path, ".fraq")) {
            throw std::runtime_error("FraqfFileWriter: expected .fraq path");
        }
        if (!ofs) throw std::runtime_error("FraqfFileWriter: open failed: " + path);

        auto fh = generate_file_header();
        ofs.write(reinterpret_cast<char*>(fh.data()), fh.size());
        if (!ofs) throw std::runtime_error("FraqfFileWriter: header write failed");
    }

    void write_compressed_block(CompressedReadBlock&& block) override {
        std::vector<char> header = generate_block_header(block);
        ofs.write(header.data(), header.size());
        ofs.write(block.name_prefix.data(), block.name_prefix.size());
        ofs.write(block.compressed_name_lengths.data(), block.compressed_name_lengths.size());
        ofs.write(block.compressed_names.data(), block.compressed_names.size());
        ofs.write(block.compressed_seq_lengths.data(), block.compressed_seq_lengths.size());
        ofs.write(block.compressed_seqs.data(), block.compressed_seqs.size());
        ofs.write(block.compressed_quals.data(), block.compressed_quals.size());
        if (!ofs) throw std::runtime_error("FraqfFileWriter: block write failed");
    }

    void flush() override {
        ofs.flush();
        if (!ofs) throw std::runtime_error("FraqfFileWriter: flush failed");
    }
    void close() override {
        ofs.close();
        if (!ofs) throw std::runtime_error("FraqfFileWriter: close failed");
    }
};

// .mem in-memory writer
struct FraqfMemWriter : IFraqfWriter {
    std::vector<CompressedReadBlock>* vec; // cached pointer to mapped vector

    FraqfMemWriter(
        const std::string& mem_key,
        std::unordered_map<std::string, std::vector<CompressedReadBlock>>& store
    ) : vec(nullptr)
    {
        if (!ends_with(mem_key, ".mem")) {
            throw std::runtime_error("FraqfMemWriter: key must end with .mem");
        }
        auto it = store.insert({mem_key, {}}).first; // create if missing
        vec = &it->second;
    }

    void write_compressed_block(CompressedReadBlock&& block) override {
        vec->push_back(std::move(block)); // destructive move into store
    }

    void flush() override {
        // nothing to do for memory sink
    }
    void close() override {
        // nothing to do for memory sink
    }
};

// Dispatcher that selects .fraq or .mem
struct FraqfMultiWriter {
    std::unique_ptr<IFraqfWriter> impl;

    // .fraq path
    explicit FraqfMultiWriter(const std::string& path) {
        if (!ends_with(path, ".fraq")) {
            throw std::runtime_error("FraqfMultiWriter: expected .fraq path");
        }
        impl = std::make_unique<FraqfFileWriter>(path);
    }

    // .mem key + store
    FraqfMultiWriter(
        const std::string& mem_key,
        std::unordered_map<std::string, std::vector<CompressedReadBlock>>& store
    ) {
        if (!ends_with(mem_key, ".mem")) {
            throw std::runtime_error("FraqfMultiWriter: expected .mem key");
        }
        impl = std::make_unique<FraqfMemWriter>(mem_key, store);
    }

    void write_compressed_block(CompressedReadBlock&& block) {
        impl->write_compressed_block(std::move(block));
    }
    void write_compressed_block(const CompressedReadBlock&) = delete; // enforce std::move

    void flush() {
        impl->flush();
    }
    void close() {
        impl->close();
    }
};


} // namespace fraq_internal

#endif // __FRAQ_WRITERS_H
