// xlsx表格读写操作
#include "utils/XlsxReader.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<uint8_t> readFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
        throw std::runtime_error("cannot open file: " + path);
    f.seekg(0, std::ios::end);
    const std::streamoff size = f.tellg();
    f.seekg(0, std::ios::beg);
    std::vector<uint8_t> data(static_cast<std::size_t>(size));
    if (size > 0)
        f.read(reinterpret_cast<char*>(data.data()), size);
    return data;
}

uint16_t le16(const uint8_t* p) {
    return static_cast<uint16_t>(p[0] | (p[1] << 8));
}

uint32_t le32(const uint8_t* p) {
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

// ---------------------- DEFLATE (inflate) ----------------------
// 参考 zlib 的 puff.c 实现，支持 stored / fixed / dynamic 三种块

const int kMaxBits = 15;
const int kMaxLcodes = 286;
const int kMaxDcodes = 30;
const int kMaxCodes = kMaxLcodes + kMaxDcodes;
const int kFixLcodes = 288;

struct Huffman {
    short count[kMaxBits + 1];
    short symbol[kMaxCodes];
};

struct Inflater {
    const uint8_t* in;
    std::size_t inlen;
    std::size_t incnt;
    uint32_t bitbuf;
    int bitcnt;
    std::vector<uint8_t>* out;

    int bits(int need) {
        uint32_t val = bitbuf;
        while (bitcnt < need) {
            if (incnt == inlen)
                throw std::runtime_error("inflate: out of input");
            val |= static_cast<uint32_t>(in[incnt++]) << bitcnt;
            bitcnt += 8;
        }
        bitbuf = val >> need;
        bitcnt -= need;
        return static_cast<int>(val & ((1u << need) - 1));
    }
};

int decode(Inflater& s, const Huffman& h) {
    int len = 1, code = 0, first = 0, index = 0;
    for (; len <= kMaxBits; ++len) {
        code |= s.bits(1);
        const int count = h.count[len];
        if (code - count < first)
            return h.symbol[index + (code - first)];
        index += count;
        first += count;
        first <<= 1;
        code <<= 1;
    }
    return -10;
}

int construct(Huffman& h, const short* length, int n) {
    for (int len = 0; len <= kMaxBits; ++len)
        h.count[len] = 0;
    for (int symbol = 0; symbol < n; ++symbol)
        h.count[length[symbol]]++;
    if (h.count[0] == n)
        return 0;
    int left = 1;
    for (int len = 1; len <= kMaxBits; ++len) {
        left <<= 1;
        left -= h.count[len];
        if (left < 0)
            return left;
    }
    short offs[kMaxBits + 1];
    offs[1] = 0;
    for (int len = 1; len < kMaxBits; ++len)
        offs[len + 1] = static_cast<short>(offs[len] + h.count[len]);
    for (int symbol = 0; symbol < n; ++symbol)
        if (length[symbol] != 0)
            h.symbol[offs[length[symbol]]++] = static_cast<short>(symbol);
    return left;
}

const short kLens[29] = {3,  4,  5,  6,  7,  8,  9,  10, 11,  13,  15,  17,
                         19, 23, 27, 31, 35, 43, 51, 59, 67,  83,  99,  115,
                         131, 163, 195, 227, 258};
const short kLext[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
                         2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
const short kDists[30] = {1,    2,    3,    4,    5,    7,    9,    13,
                          17,   25,   33,   49,   65,   97,   129,  193,
                          257,  385,  513,  769,  1025, 1537, 2049, 3073,
                          4097, 6145, 8193, 12289, 16385, 24577};
const short kDext[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6,
                         6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

int codes(Inflater& s, const Huffman& lencode, const Huffman& distcode) {
    for (;;) {
        int symbol = decode(s, lencode);
        if (symbol < 0)
            return symbol;
        if (symbol < 256) {
            s.out->push_back(static_cast<uint8_t>(symbol));
        } else if (symbol == 256) {
            return 0;
        } else {
            symbol -= 257;
            if (symbol >= 29)
                return -9;
            const int len = kLens[symbol] + s.bits(kLext[symbol]);
            const int dsym = decode(s, distcode);
            if (dsym < 0)
                return dsym;
            const int dist = kDists[dsym] + s.bits(kDext[dsym]);
            if (static_cast<std::size_t>(dist) > s.out->size())
                return -11;
            for (int i = 0; i < len; ++i)
                s.out->push_back((*s.out)[s.out->size() - dist]);
        }
    }
}

int stored(Inflater& s) {
    s.bitbuf = 0;
    s.bitcnt = 0;
    if (s.incnt + 4 > s.inlen)
        return 2;
    const int len = s.in[s.incnt] | (s.in[s.incnt + 1] << 8);
    const int nlen = s.in[s.incnt + 2] | (s.in[s.incnt + 3] << 8);
    s.incnt += 4;
    if ((len ^ 0xffff) != nlen)
        return -4;
    if (s.incnt + len > s.inlen)
        return 2;
    for (int i = 0; i < len; ++i)
        s.out->push_back(s.in[s.incnt++]);
    return 0;
}

int fixed(Inflater& s) {
    Huffman lencode, distcode;
    short lengths[kFixLcodes];
    int symbol;
    for (symbol = 0; symbol < 144; ++symbol)
        lengths[symbol] = 8;
    for (; symbol < 256; ++symbol)
        lengths[symbol] = 9;
    for (; symbol < 280; ++symbol)
        lengths[symbol] = 7;
    for (; symbol < kFixLcodes; ++symbol)
        lengths[symbol] = 8;
    // 固定表是合法的：construct 对不完整码表(定长距离码)会返回正的 left，
    // 这里只需在负值(过订)时才报错，其余情况忽略返回值。
    int err = construct(lencode, lengths, kFixLcodes);
    if (err < 0)
        return err;
    for (symbol = 0; symbol < kMaxDcodes; ++symbol)
        lengths[symbol] = 5;
    err = construct(distcode, lengths, kMaxDcodes);
    if (err < 0)
        return err;
    return codes(s, lencode, distcode);
}

int dynamic(Inflater& s) {
    static const short kOrder[19] = {16, 17, 18, 0, 8,  7, 9,  6, 10, 5,
                                     11, 4,  12, 3, 13, 2, 14, 1, 15};
    Huffman lencode, distcode;
    short lengths[kMaxCodes];
    for (int i = 0; i < kMaxCodes; ++i)
        lengths[i] = 0;

    const int nlen = s.bits(5) + 257;
    const int ndist = s.bits(5) + 1;
    const int ncode = s.bits(4) + 4;
    if (nlen > kMaxLcodes || ndist > kMaxDcodes)
        return -3;

    short lens[19];
    for (int i = 0; i < 19; ++i)
        lens[i] = 0;
    int index = 0;
    for (; index < ncode; ++index)
        lens[kOrder[index]] = static_cast<short>(s.bits(3));

    int err = construct(lencode, lens, 19);
    if (err)
        return -4;

    index = 0;
    while (index < nlen + ndist) {
        const int symbol = decode(s, lencode);
        if (symbol < 0)
            return symbol;
        if (symbol < 16) {
            lengths[index++] = static_cast<short>(symbol);
        } else {
            int len = 0;
            int count;
            if (symbol == 16) {
                if (index == 0)
                    return -5;
                len = lengths[index - 1];
                count = 3 + s.bits(2);
            } else if (symbol == 17) {
                count = 3 + s.bits(3);
            } else {
                count = 11 + s.bits(7);
            }
            if (index + count > nlen + ndist)
                return -6;
            while (count--)
                lengths[index++] = static_cast<short>(len);
        }
    }
    if (lengths[256] == 0)
        return -7;

    err = construct(lencode, lengths, nlen);
    if (err && (err < 0 || nlen != lencode.count[0] + lencode.count[1]))
        return -8;
    err = construct(distcode, lengths + nlen, ndist);
    if (err && (err < 0 || ndist != distcode.count[0] + distcode.count[1]))
        return -9;
    return codes(s, lencode, distcode);
}

std::vector<uint8_t> inflate(const std::vector<uint8_t>& in,
                             std::size_t uncompSizeHint = 0) {
    Inflater s{in.data(), in.size(), 0, 0, 0, nullptr};
    std::vector<uint8_t> out;
    // 按中央目录记录的未压缩大小预分配，避免逐字节扩容
    if (uncompSizeHint > 0)
        out.reserve(uncompSizeHint);
    s.out = &out;
    int last;
    do {
        last = s.bits(1);
        const int type = s.bits(2);
        const int err = (type == 0)   ? stored(s)
                        : (type == 1) ? fixed(s)
                        : (type == 2) ? dynamic(s)
                                      : -1;
        if (err != 0)
            throw std::runtime_error("inflate failed");
    } while (!last);
    return out;
}

// ---------------------- DEFLATE (deflate) ----------------------
// 仅使用固定 Huffman 码 + 贪心 LZ77 匹配，输出原始(无 zlib 头) deflate 流。
// 目标是把 xlsx 内文本(XML)压得比 stored 小，且实现足够简单可靠。

struct BitWriter {
    std::vector<uint8_t>& out;
    uint32_t bitbuf;
    int bitcnt;

    explicit BitWriter(std::vector<uint8_t>& o)
        : out(o), bitbuf(0), bitcnt(0) {}

    // 写入 nbits 位原始整数(LSB-first)，用于额外位
    void put(uint32_t value, int nbits) {
        bitbuf |= (value & ((1u << nbits) - 1u)) << bitcnt;
        bitcnt += nbits;
        while (bitcnt >= 8) {
            out.push_back(static_cast<uint8_t>(bitbuf & 0xFF));
            bitbuf >>= 8;
            bitcnt -= 8;
        }
    }

    // 写入 Huffman 码(MSB-first)
    void putHuff(uint32_t code, int nbits) {
        for (int i = nbits - 1; i >= 0; --i)
            put((code >> i) & 1u, 1);
    }

    void flush() {
        if (bitcnt > 0) {
            out.push_back(static_cast<uint8_t>(bitbuf & 0xFF));
            bitbuf = 0;
            bitcnt = 0;
        }
    }
};

// 固定 Huffman 表中的 literal/length 码
void emitLitLen(BitWriter& w, int sym) {
    if (sym <= 143)
        w.putHuff(static_cast<uint32_t>(0x30 + sym), 8);
    else if (sym <= 255)
        w.putHuff(static_cast<uint32_t>(0x190 + (sym - 144)), 9);
    else if (sym <= 279)
        w.putHuff(static_cast<uint32_t>(sym - 256), 7);
    else
        w.putHuff(static_cast<uint32_t>(0xC0 + (sym - 280)), 8);
}

void emitMatch(BitWriter& w, int len, int dist) {
    int lsym = 0;
    for (int i = 28; i >= 0; --i)
        if (len >= kLens[i]) {
            lsym = i;
            break;
        }
    emitLitLen(w, 257 + lsym);
    w.put(static_cast<uint32_t>(len - kLens[lsym]), kLext[lsym]);

    int dsym = 0;
    for (int i = 29; i >= 0; --i)
        if (dist >= kDists[i]) {
            dsym = i;
            break;
        }
    w.putHuff(static_cast<uint32_t>(dsym), 5);
    w.put(static_cast<uint32_t>(dist - kDists[dsym]), kDext[dsym]);
}

std::vector<uint8_t> deflate(const std::vector<uint8_t>& in) {
    std::vector<uint8_t> out;
    out.reserve(in.size() / 2 + 16);
    BitWriter w(out);
    w.put(1, 1);  // BFINAL = 1
    w.put(1, 2);  // BTYPE  = 01 (fixed Huffman)

    const size_t n = in.size();
    const int kHashBits = 15;
    const int kHashSize = 1 << kHashBits;
    std::vector<int> head(static_cast<std::size_t>(kHashSize), -1);
    std::vector<int> prev(n, -1);

    auto hash3 = [&in](std::size_t p) -> int {
        uint32_t h = static_cast<uint32_t>(in[p]) |
                     (static_cast<uint32_t>(in[p + 1]) << 8) |
                     (static_cast<uint32_t>(in[p + 2]) << 16);
        h *= 2654435761u;
        return static_cast<int>(h >> (32 - kHashBits));
    };

    std::size_t i = 0;
    while (i < n) {
        int bestLen = 0;
        int bestDist = 0;
        if (i + 3 <= n) {
            const int h = hash3(i);
            int cand = head[h];
            int chain = 0;
            const int kMaxChain = 128;
            std::size_t maxLen = n - i;
            if (maxLen > 258)
                maxLen = 258;
            while (cand >= 0 && chain++ < kMaxChain &&
                   static_cast<int>(i) - cand <= 32768) {
                std::size_t l = 0;
                while (l < maxLen && in[static_cast<std::size_t>(cand) + l] ==
                                          in[i + l])
                    ++l;
                if (static_cast<int>(l) > bestLen) {
                    bestLen = static_cast<int>(l);
                    bestDist = static_cast<int>(i) - cand;
                    if (l >= maxLen)
                        break;
                }
                cand = prev[static_cast<std::size_t>(cand)];
            }
            prev[i] = head[h];
            head[h] = static_cast<int>(i);
        }

        if (bestLen >= 3) {
            emitMatch(w, bestLen, bestDist);
            for (int k = 1; k < bestLen; ++k) {
                const std::size_t p = i + static_cast<std::size_t>(k);
                if (p + 3 <= n) {
                    const int hh = hash3(p);
                    prev[p] = head[hh];
                    head[hh] = static_cast<int>(p);
                }
            }
            i += static_cast<std::size_t>(bestLen);
        } else {
            emitLitLen(w, in[i]);
            ++i;
        }
    }
    emitLitLen(w, 256);  // end of block
    w.flush();
    return out;
}

// ---------------------- ZIP 解析 ----------------------

struct ZipEntry {
    std::string name;
    uint16_t method;
    uint32_t crc;
    uint32_t compSize;    // 压缩后大小(数据区字节数)
    uint32_t uncompSize;  // 解压后大小
    uint32_t localOffset;
};

std::vector<ZipEntry> readCentralDirectory(const std::vector<uint8_t>& d) {
    if (d.size() < 22)
        throw std::runtime_error("file too small to be xlsx");
    std::size_t eocd = std::string::npos;
    const std::size_t lo = d.size() >= 65557 ? d.size() - 65557 : 0;
    for (std::size_t i = d.size() - 22 + 1; i-- > lo;) {
        if (le32(&d[i]) == 0x06054b50u) {
            eocd = i;
            break;
        }
    }
    if (eocd == std::string::npos)
        throw std::runtime_error("end of central directory not found");

    const uint32_t count = le16(&d[eocd + 10]);
    const uint32_t cdOffset = le32(&d[eocd + 16]);

    std::vector<ZipEntry> entries;
    std::size_t p = cdOffset;
    for (uint32_t k = 0; k < count; ++k) {
        if (p + 46 > d.size() || le32(&d[p]) != 0x02014b50u)
            break;
        ZipEntry e;
        e.method = le16(&d[p + 10]);
        e.crc = le32(&d[p + 16]);
        e.compSize = le32(&d[p + 20]);
        e.uncompSize = le32(&d[p + 24]);
        const uint16_t nameLen = le16(&d[p + 28]);
        const uint16_t extraLen = le16(&d[p + 30]);
        const uint16_t commentLen = le16(&d[p + 32]);
        e.localOffset = le32(&d[p + 42]);
        e.name.assign(reinterpret_cast<const char*>(&d[p + 46]), nameLen);
        entries.push_back(e);
        p += 46 + nameLen + extraLen + commentLen;
    }
    return entries;
}

std::vector<uint8_t> extractEntry(const std::vector<uint8_t>& d,
                                  const ZipEntry& e) {
    const std::size_t p = e.localOffset;
    if (p + 30 > d.size() || le32(&d[p]) != 0x04034b50u)
        throw std::runtime_error("bad local file header");
    const uint16_t nameLen = le16(&d[p + 26]);
    const uint16_t extraLen = le16(&d[p + 28]);
    const std::size_t dataOff = p + 30 + nameLen + extraLen;
    if (dataOff + e.compSize > d.size())
        throw std::runtime_error("zip data out of range");
    std::vector<uint8_t> raw(d.begin() + dataOff,
                             d.begin() + dataOff + e.compSize);
    if (e.method == 0)
        return raw;
    if (e.method == 8)
        return inflate(raw, e.uncompSize);
    throw std::runtime_error("unsupported zip compression method");
}

const ZipEntry* findEntry(const std::vector<ZipEntry>& entries,
                          const std::string& name) {
    for (const auto& e : entries)
        if (e.name == name)
            return &e;
    return nullptr;
}

std::string toString(const std::vector<uint8_t>& bytes) {
    return std::string(bytes.begin(), bytes.end());
}

// ---------------------- XML 工具 ----------------------

std::string getAttr(const std::string& element, const std::string& name) {
    const std::string key = name + "=\"";
    const std::size_t p = element.find(key);
    if (p == std::string::npos)
        return "";
    const std::size_t begin = p + key.size();
    const std::size_t q = element.find('"', begin);
    if (q == std::string::npos)
        return "";
    return element.substr(begin, q - begin);
}

std::string unescapeXml(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '&') {
            if (s.compare(i, 5, "&amp;") == 0) {
                out += '&';
                i += 4;
            } else if (s.compare(i, 4, "&lt;") == 0) {
                out += '<';
                i += 3;
            } else if (s.compare(i, 4, "&gt;") == 0) {
                out += '>';
                i += 3;
            } else if (s.compare(i, 6, "&quot;") == 0) {
                out += '"';
                i += 5;
            } else if (s.compare(i, 6, "&apos;") == 0) {
                out += '\'';
                i += 5;
            } else {
                out += s[i];
            }
        } else {
            out += s[i];
        }
    }
    return out;
}

// 拼接元素内所有 <t>...</t> 的文本(用于共享字符串 / 内联字符串)
std::string extractText(const std::string& s) {
    std::string out;
    std::size_t p = 0;
    for (;;) {
        const std::size_t t = s.find("<t", p);
        if (t == std::string::npos)
            break;
        const std::size_t gt = s.find('>', t);
        if (gt == std::string::npos)
            break;
        const std::size_t close = s.find("</t>", gt);
        if (close == std::string::npos)
            break;
        out += s.substr(gt + 1, close - (gt + 1));
        p = close + 4;
    }
    return unescapeXml(out);
}

// 取 <v>...</v> 的原始内容
std::string extractRawV(const std::string& cell) {
    const std::size_t vp = cell.find("<v>");
    if (vp == std::string::npos)
        return "";
    const std::size_t vq = cell.find("</v>", vp);
    if (vq == std::string::npos)
        return "";
    return cell.substr(vp + 3, vq - (vp + 3));
}

// ---------------------- 工作表解析 ----------------------

struct Cell {
    bool hasValue = false;
    bool isString = false;
    double number = 0.0;
    std::string text;
};

// 解析共享字符串表 sharedStrings.xml，返回按索引排列的字符串
std::vector<std::string> parseSharedStrings(const std::string& xml) {
    std::vector<std::string> shared;
    std::size_t p = 0;
    for (;;) {
        const std::size_t s = xml.find("<si>", p);
        if (s == std::string::npos)
            break;
        const std::size_t e = xml.find("</si>", s);
        if (e == std::string::npos)
            break;
        shared.push_back(extractText(xml.substr(s, e - s)));
        p = e + 5;
    }
    return shared;
}

// 定位下一个单元格起始(<c ...> / <c> / <c/>)，
// 只认标签名恰为 "c" 的位置，避免命中 <conditionalFormatting 等以 c 开头的标签。
std::size_t findCellStart(const std::string& s, std::size_t from) {
    std::size_t p = from;
    while ((p = s.find("<c", p)) != std::string::npos) {
        const char n = (p + 2 < s.size()) ? s[p + 2] : '\0';
        if (n == ' ' || n == '>' || n == '/' || n == '\t' || n == '\n' ||
            n == '\r')
            return p;
        p += 2;
    }
    return std::string::npos;
}

// 从引用串(如 B3 / $B$3)解析行列号，行列均 1 起始
bool parseRef(const std::string& ref, int& row, int& col) {
    row = 0;
    col = 0;
    for (const char ch : ref) {
        if (ch == '$')
            continue;
        if (ch >= 'A' && ch <= 'Z')
            col = col * 26 + (ch - 'A' + 1);
        else if (ch >= 'a' && ch <= 'z')
            col = col * 26 + (ch - 'a' + 1);
        else if (ch >= '0' && ch <= '9')
            row = row * 10 + (ch - '0');
    }
    return col > 0 && row > 0;
}

// 把一个单元格 XML 解码为 Cell(内联串/共享串/公式串/布尔/错误/数值)
Cell decodeCell(const std::string& cell, const std::vector<std::string>& shared) {
    Cell value;
    if (cell.find("t=\"inlineStr\"") != std::string::npos) {
        value.hasValue = true;
        value.isString = true;
        value.text = extractText(cell);
    } else if (cell.find("t=\"s\"") != std::string::npos) {
        const std::string v = extractRawV(cell);
        if (!v.empty()) {
            try {
                const std::size_t idx =
                    static_cast<std::size_t>(std::stoul(v));
                if (idx < shared.size()) {
                    value.hasValue = true;
                    value.isString = true;
                    value.text = shared[idx];
                }
            } catch (...) {
            }
        }
    } else if (cell.find("t=\"str\"") != std::string::npos) {
        value.hasValue = true;
        value.isString = true;
        value.text = unescapeXml(extractRawV(cell));
    } else if (cell.find("t=\"b\"") != std::string::npos) {
        const std::string v = extractRawV(cell);
        if (!v.empty()) {
            value.hasValue = true;
            value.isString = false;
            value.number = (v == "1") ? 1.0 : 0.0;
        }
    } else if (cell.find("t=\"e\"") != std::string::npos) {
        // 错误值，忽略
    } else {
        const std::string v = extractRawV(cell);
        if (!v.empty()) {
            try {
                value.number = std::stod(v);
                value.hasValue = true;
            } catch (...) {
            }
        }
    }
    return value;
}

// 解析工作表: 返回 cells[行][列] = Cell，行列均为 1 起始
std::map<int, std::map<int, Cell>> parseSheet(
    const std::string& xml, const std::vector<std::string>& shared) {
    std::map<int, std::map<int, Cell>> cells;
    std::size_t pos = 0;
    for (;;) {
        const std::size_t c = findCellStart(xml, pos);
        if (c == std::string::npos)
            break;
        // 先定位 <c ...> 起始标签的结束位置，据其是否以 "/" 结尾判断单元格是否自闭合，
        // 避免把单元格内部的 <f .../> 误当作单元格结束而丢掉后面的 <v>
        const std::size_t tagEnd = xml.find('>', c);
        if (tagEnd == std::string::npos)
            break;
        std::size_t end;
        std::size_t next;
        if (tagEnd > c && xml[tagEnd - 1] == '/') {
            end = tagEnd - 1;
            next = tagEnd + 1;
        } else {
            const std::size_t endClose = xml.find("</c>", tagEnd);
            if (endClose == std::string::npos)
                break;
            end = endClose;
            next = endClose + 4;
        }
        const std::string cell = xml.substr(c, end - c);
        pos = next;

        const std::string ref = getAttr(cell, "r");
        if (ref.empty())
            continue;

        int row = 0;
        int col = 0;
        if (!parseRef(ref, row, col))
            continue;

        Cell value = decodeCell(cell, shared);
        if (value.hasValue)
            cells[row][col] = value;
    }
    return cells;
}

// 由 workbook.xml 与 workbook.xml.rels 求工作表对应的 zip 路径
std::string resolveSheetPath(const std::vector<uint8_t>& data,
                             const std::vector<ZipEntry>& entries,
                             const std::string& sheetName) {
    const ZipEntry* wb = findEntry(entries, "xl/workbook.xml");
    if (wb == nullptr)
        throw std::runtime_error("xl/workbook.xml not found");
    const std::string wbxml = toString(extractEntry(data, *wb));

    std::string rid;
    std::size_t p = 0;
    for (;;) {
        const std::size_t s = wbxml.find("<sheet ", p);
        if (s == std::string::npos)
            break;
        const std::size_t e = wbxml.find('>', s);
        if (e == std::string::npos)
            break;
        const std::string el = wbxml.substr(s, e - s);
        if (getAttr(el, "name") == sheetName) {
            rid = getAttr(el, "r:id");
            break;
        }
        p = e + 1;
    }
    if (rid.empty())
        throw std::runtime_error("sheet not found: " + sheetName);

    const ZipEntry* rel = findEntry(entries, "xl/_rels/workbook.xml.rels");
    if (rel == nullptr)
        throw std::runtime_error("xl/_rels/workbook.xml.rels not found");
    const std::string relxml = toString(extractEntry(data, *rel));

    std::string target;
    p = 0;
    for (;;) {
        const std::size_t s = relxml.find("<Relationship ", p);
        if (s == std::string::npos)
            break;
        const std::size_t e = relxml.find('>', s);
        if (e == std::string::npos)
            break;
        const std::string el = relxml.substr(s, e - s);
        if (getAttr(el, "Id") == rid) {
            target = getAttr(el, "Target");
            break;
        }
        p = e + 1;
    }
    if (target.empty())
        throw std::runtime_error("relationship not found for sheet: " +
                                 sheetName);

    if (target[0] == '/')
        return target.substr(1);
    return "xl/" + target;
}

// ---------------------- 写出 (ZIP / worksheet) ----------------------

// 一个待写出的 zip 部件。未修改的部件保留其原始压缩字节与 method，
// 只有被重建的部件才重新压缩/存储。
struct OutPart {
    std::string name;
    uint16_t method = 0;        // 0 stored, 8 deflate
    uint32_t crc = 0;           // 未压缩数据的 crc32
    uint32_t compSize = 0;      // data 中的字节数
    uint32_t uncompSize = 0;
    std::vector<uint8_t> data;  // 实际写入 zip 的字节(已压缩或原样)
};

uint32_t crc32(const uint8_t* data, std::size_t len) {
    uint32_t crc = 0xFFFFFFFFu;
    for (std::size_t i = 0; i < len; ++i) {
        crc ^= data[i];
        for (int j = 0; j < 8; ++j)
            crc = (crc & 1u) ? ((crc >> 1) ^ 0xEDB88320u) : (crc >> 1);
    }
    return ~crc;
}

std::string colName(int col) {
    std::string s;
    while (col > 0) {
        const int rem = (col - 1) % 26;
        s.insert(s.begin(), static_cast<char>('A' + rem));
        col = (col - 1) / 26;
    }
    return s;
}

std::string escapeXml(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (const char ch : s) {
        switch (ch) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out += ch;
        }
    }
    return out;
}

std::string formatNumber(double v) {
    std::ostringstream os;
    os << std::setprecision(15) << v;
    return os.str();
}

void appendLe16(std::vector<uint8_t>& b, uint16_t v) {
    b.push_back(static_cast<uint8_t>(v & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 8) & 0xff));
}

void appendLe32(std::vector<uint8_t>& b, uint32_t v) {
    b.push_back(static_cast<uint8_t>(v & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 8) & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 16) & 0xff));
    b.push_back(static_cast<uint8_t>((v >> 24) & 0xff));
}

std::vector<uint8_t> toBytes(const std::string& s) {
    return std::vector<uint8_t>(s.begin(), s.end());
}

std::string bytesToString(const std::vector<uint8_t>& b) {
    return std::string(b.begin(), b.end());
}

// 从完整 zip 数据中取出某条目未经解压的原始字节(用于原样保留)
std::vector<uint8_t> rawEntryBytes(const std::vector<uint8_t>& d,
                                   const ZipEntry& e) {
    const std::size_t p = e.localOffset;
    if (p + 30 > d.size() || le32(&d[p]) != 0x04034b50u)
        throw std::runtime_error("bad local file header");
    const uint16_t nameLen = le16(&d[p + 26]);
    const uint16_t extraLen = le16(&d[p + 28]);
    const std::size_t dataOff = p + 30 + nameLen + extraLen;
    if (dataOff + e.compSize > d.size())
        throw std::runtime_error("zip data out of range");
    return std::vector<uint8_t>(d.begin() + dataOff,
                                d.begin() + dataOff + e.compSize);
}

// 压缩内容；若压缩后反而更大则退回 stored
OutPart makeDeflated(const std::string& name, const std::string& content) {
    const std::vector<uint8_t> raw = toBytes(content);
    OutPart p;
    p.name = name;
    p.uncompSize = static_cast<uint32_t>(raw.size());
    p.crc = crc32(raw.data(), raw.size());
    const std::vector<uint8_t> comp = deflate(raw);
    if (comp.size() < raw.size()) {
        p.method = 8;
        p.data = comp;
    } else {
        p.method = 0;
        p.data = raw;
    }
    p.compSize = static_cast<uint32_t>(p.data.size());
    return p;
}

// 取部件未压缩内容
std::string partText(const OutPart& p) {
    if (p.method == 0)
        return bytesToString(p.data);
    if (p.method == 8)
        return bytesToString(inflate(p.data, p.uncompSize));
    throw std::runtime_error("unsupported method in part");
}

// 按各部件自身的压缩方式写出 zip
void writeZip(const std::string& path, const std::vector<OutPart>& parts) {
    struct Central {
        std::string name;
        uint16_t method;
        uint32_t crc;
        uint32_t comp;
        uint32_t uncomp;
        uint32_t offset;
    };
    std::vector<uint8_t> out;
    std::vector<Central> central;

    for (const auto& part : parts) {
        Central c{part.name, part.method, part.crc, part.compSize,
                  part.uncompSize, static_cast<uint32_t>(out.size())};

        appendLe32(out, 0x04034b50u);
        appendLe16(out, 20);
        appendLe16(out, 0);
        appendLe16(out, part.method);
        appendLe16(out, 0);     // time
        appendLe16(out, 0x21);  // date (1980-01-01)
        appendLe32(out, part.crc);
        appendLe32(out, part.compSize);
        appendLe32(out, part.uncompSize);
        appendLe16(out, static_cast<uint16_t>(part.name.size()));
        appendLe16(out, 0);
        out.insert(out.end(), part.name.begin(), part.name.end());
        out.insert(out.end(), part.data.begin(), part.data.end());
        central.push_back(c);
    }

    const uint32_t cdOffset = static_cast<uint32_t>(out.size());
    for (const auto& c : central) {
        appendLe32(out, 0x02014b50u);
        appendLe16(out, 20);
        appendLe16(out, 20);
        appendLe16(out, 0);
        appendLe16(out, c.method);
        appendLe16(out, 0);
        appendLe16(out, 0x21);
        appendLe32(out, c.crc);
        appendLe32(out, c.comp);
        appendLe32(out, c.uncomp);
        appendLe16(out, static_cast<uint16_t>(c.name.size()));
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe16(out, 0);
        appendLe32(out, 0);
        appendLe32(out, c.offset);
        out.insert(out.end(), c.name.begin(), c.name.end());
    }
    const uint32_t cdSize = static_cast<uint32_t>(out.size()) - cdOffset;

    appendLe32(out, 0x06054b50u);
    appendLe16(out, 0);
    appendLe16(out, 0);
    appendLe16(out, static_cast<uint16_t>(central.size()));
    appendLe16(out, static_cast<uint16_t>(central.size()));
    appendLe32(out, cdSize);
    appendLe32(out, cdOffset);
    appendLe16(out, 0);

    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    if (!f)
        throw std::runtime_error("cannot write file: " + path);
    f.write(reinterpret_cast<const char*>(out.data()),
            static_cast<std::streamsize>(out.size()));
}

std::vector<OutPart> buildEmptyPackage() {
    std::vector<OutPart> parts;
    parts.push_back(makeDeflated(
        "[Content_Types].xml",
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/"
        "content-types\"><Default Extension=\"rels\" "
        "ContentType=\"application/vnd.openxmlformats-package."
        "relationships+xml\"/><Default Extension=\"xml\" "
        "ContentType=\"application/xml\"/><Override "
        "PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd."
        "openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
        "</Types>"));
    parts.push_back(makeDeflated(
        "_rels/.rels",
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/"
        "package/2006/relationships\"><Relationship Id=\"rId1\" "
        "Type=\"http://schemas.openxmlformats.org/officeDocument/2006/"
        "relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
        "</Relationships>"));
    parts.push_back(makeDeflated(
        "xl/workbook.xml",
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<workbook xmlns=\"http://schemas.openxmlformats.org/"
        "spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats."
        "org/officeDocument/2006/relationships\"><sheets></sheets>"
        "</workbook>"));
    parts.push_back(makeDeflated(
        "xl/_rels/workbook.xml.rels",
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/"
        "package/2006/relationships\"></Relationships>"));
    return parts;
}

void replacePart(std::vector<OutPart>& parts, const std::string& name,
                 const OutPart& np) {
    for (auto& p : parts)
        if (p.name == name) {
            p = np;
            return;
        }
    parts.push_back(np);
}

// ---------------------- 工作表增量修补 ----------------------

struct RawCell {
    std::string raw;  // 原始 <c ...>...</c> 或 <c .../>
    Cell value;
};

struct RawRow {
    std::string openTag;  // 原始 <row ...>；合成行为空串
    bool selfClosed = false;
    std::map<int, RawCell> cells;
};

struct RawSheet {
    bool hasSheetData = false;
    std::string prefix;  // <sheetData ...> 之前
    std::string suffix;  // </sheetData> 之后
    std::map<int, RawRow> rows;
};

RawSheet parseRawSheet(const std::string& xml,
                       const std::vector<std::string>& shared) {
    RawSheet rs;
    const std::size_t sd = xml.find("<sheetData");
    if (sd == std::string::npos)
        return rs;
    const std::size_t sdEnd = xml.find('>', sd);
    if (sdEnd == std::string::npos)
        return rs;
    rs.hasSheetData = true;
    rs.prefix = xml.substr(0, sd);

    const bool sdSelfClosed = (xml[sdEnd - 1] == '/');
    std::string inner;
    if (sdSelfClosed) {
        rs.suffix = xml.substr(sdEnd + 1);
    } else {
        const std::size_t close = xml.find("</sheetData>", sdEnd);
        if (close == std::string::npos)
            return rs;
        inner = xml.substr(sdEnd + 1, close - (sdEnd + 1));
        rs.suffix = xml.substr(close + 12);  // strlen("</sheetData>")
    }

    std::size_t p = 0;
    int autoRow = 0;
    while (true) {
        const std::size_t r = inner.find("<row", p);
        if (r == std::string::npos)
            break;
        const char nc = (r + 4 < inner.size()) ? inner[r + 4] : '\0';
        if (!(nc == ' ' || nc == '>' || nc == '/' || nc == '\t' ||
              nc == '\n' || nc == '\r')) {
            p = r + 4;
            continue;
        }
        const std::size_t tagEnd = inner.find('>', r);
        if (tagEnd == std::string::npos)
            break;
        const bool rowSelf = (inner[tagEnd - 1] == '/');
        std::size_t next;
        std::string rowInner;
        if (rowSelf) {
            next = tagEnd + 1;
        } else {
            const std::size_t rowEnd = inner.find("</row>", tagEnd);
            if (rowEnd == std::string::npos)
                break;
            rowInner = inner.substr(tagEnd + 1, rowEnd - (tagEnd + 1));
            next = rowEnd + 6;
        }
        const std::string openTag = inner.substr(r, tagEnd - r + 1);

        int rowNum = autoRow + 1;
        const std::string num = getAttr(openTag, "r");
        if (!num.empty()) {
            try {
                rowNum = std::stoi(num);
            } catch (...) {
            }
        }
        autoRow = rowNum;

        RawRow rr;
        rr.openTag = openTag;
        rr.selfClosed = rowSelf;
        if (!rowSelf) {
            std::size_t cp = 0;
            while (true) {
                const std::size_t c = findCellStart(rowInner, cp);
                if (c == std::string::npos)
                    break;
                const std::size_t cTagEnd = rowInner.find('>', c);
                if (cTagEnd == std::string::npos)
                    break;
                std::size_t cEnd, cNext;
                if (cTagEnd > c && rowInner[cTagEnd - 1] == '/') {
                    cEnd = cTagEnd + 1;
                    cNext = cTagEnd + 1;
                } else {
                    const std::size_t cClose = rowInner.find("</c>", cTagEnd);
                    if (cClose == std::string::npos)
                        break;
                    cEnd = cClose + 4;
                    cNext = cEnd;
                }
                const std::string rawCell = rowInner.substr(c, cEnd - c);
                cp = cNext;

                const std::string ref = getAttr(rawCell, "r");
                int rrow = 0, ccol = 0;
                if (ref.empty() || !parseRef(ref, rrow, ccol))
                    continue;
                RawCell rc;
                rc.raw = rawCell;
                rc.value = decodeCell(rawCell, shared);
                rr.cells[ccol] = rc;
            }
        }
        rs.rows[rowNum] = rr;
        p = next;
    }
    return rs;
}

std::string serializeRows(const std::map<int, RawRow>& rows) {
    std::string out;
    for (const auto& kv : rows) {
        const int num = kv.first;
        const RawRow& r = kv.second;
        std::string cellsXml;
        for (const auto& c : r.cells)
            cellsXml += c.second.raw;
        if (r.openTag.empty()) {
            out += "<row r=\"" + std::to_string(num) + "\">" + cellsXml +
                   "</row>";
        } else if (r.selfClosed && cellsXml.empty()) {
            out += r.openTag;
        } else {
            std::string tag = r.openTag;
            if (r.selfClosed)
                tag = tag.substr(0, tag.size() - 2) + ">";
            out += tag + cellsXml + "</row>";
        }
    }
    return out;
}

// 在已有工作表 XML 上就地增补表头/数据，保留其余单元格与格式
std::string patchSheetXml(const std::string& sheetXml,
                          const std::vector<std::string>& shared,
                          const std::vector<std::string>& headers,
                          const std::vector<std::vector<double>>& columns,
                          bool overwrite, int explicitStartRow) {
    RawSheet rs = parseRawSheet(sheetXml, shared);

    // 在指定行写一个新块：表头在 explicitStartRow，数据从下一行开始，
    // 保留工作表中其它已有单元格。
    if (explicitStartRow > 0) {
        const std::size_t rows = columns.empty() ? 0 : columns[0].size();
        if (!rs.hasSheetData) {
            std::string s =
                "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
                "<worksheet xmlns=\"http://schemas.openxmlformats.org/"
                "spreadsheetml/2006/main\"><sheetData>";
            s += "<row r=\"" + std::to_string(explicitStartRow) + "\">";
            for (std::size_t i = 0; i < headers.size(); ++i)
                s += "<c r=\"" + colName(static_cast<int>(i) + 1) +
                     std::to_string(explicitStartRow) +
                     "\" t=\"inlineStr\"><is><t>" + escapeXml(headers[i]) +
                     "</t></is></c>";
            s += "</row>";
            for (std::size_t r = 0; r < rows; ++r) {
                const int rowNum = explicitStartRow + 1 + static_cast<int>(r);
                s += "<row r=\"" + std::to_string(rowNum) + "\">";
                for (std::size_t i = 0; i < headers.size(); ++i) {
                    if (std::isnan(columns[i][r]))
                        continue;
                    s += "<c r=\"" + colName(static_cast<int>(i) + 1) +
                         std::to_string(rowNum) + "\"><v>" +
                         formatNumber(columns[i][r]) + "</v></c>";
                }
                s += "</row>";
            }
            s += "</sheetData></worksheet>";
            return s;
        }

        for (std::size_t i = 0; i < headers.size(); ++i) {
            RawCell rc;
            rc.raw = "<c r=\"" + colName(static_cast<int>(i) + 1) +
                     std::to_string(explicitStartRow) +
                     "\" t=\"inlineStr\"><is><t>" + escapeXml(headers[i]) +
                     "</t></is></c>";
            rc.value.hasValue = true;
            rc.value.isString = true;
            rc.value.text = headers[i];
            rs.rows[explicitStartRow].cells[static_cast<int>(i) + 1] = rc;
        }
        for (std::size_t r = 0; r < rows; ++r) {
            const int rowNum = explicitStartRow + 1 + static_cast<int>(r);
            for (std::size_t i = 0; i < headers.size(); ++i) {
                if (std::isnan(columns[i][r]))
                    continue;
                RawCell rc;
                rc.raw = "<c r=\"" + colName(static_cast<int>(i) + 1) +
                         std::to_string(rowNum) + "\"><v>" +
                         formatNumber(columns[i][r]) + "</v></c>";
                rc.value.hasValue = true;
                rc.value.isString = false;
                rc.value.number = columns[i][r];
                rs.rows[rowNum].cells[static_cast<int>(i) + 1] = rc;
            }
        }
        return rs.prefix + "<sheetData>" + serializeRows(rs.rows) +
               "</sheetData>" + rs.suffix;
    }

    if (overwrite || !rs.hasSheetData) {
        // 无 sheetData：退化为重建一张最小工作表
        std::string s =
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
            "<worksheet xmlns=\"http://schemas.openxmlformats.org/"
            "spreadsheetml/2006/main\"><sheetData>";
        s += "<row r=\"1\">";
        for (std::size_t i = 0; i < headers.size(); ++i)
            s += "<c r=\"" + colName(static_cast<int>(i) + 1) +
                 "1\" t=\"inlineStr\"><is><t>" + escapeXml(headers[i]) +
                 "</t></is></c>";
        s += "</row>";
        const std::size_t rows = columns.empty() ? 0 : columns[0].size();
        for (std::size_t r = 0; r < rows; ++r) {
            s += "<row r=\"" + std::to_string(r + 2) + "\">";
            for (std::size_t i = 0; i < headers.size(); ++i) {
                if (std::isnan(columns[i][r]))
                    continue;
                s += "<c r=\"" + colName(static_cast<int>(i) + 1) +
                     std::to_string(r + 2) + "\"><v>" +
                     formatNumber(columns[i][r]) + "</v></c>";
            }
            s += "</row>";
        }
        s += "</sheetData></worksheet>";
        return s;
    }

    // 表头行：出现任一目标表头名的行；都没有则用第 1 行
    int headerRow = 1;
    bool found = false;
    for (const auto& rk : rs.rows) {
        for (const auto& ck : rk.second.cells) {
            const Cell& cell = ck.second.value;
            if (!cell.isString)
                continue;
            for (const std::string& h : headers)
                if (cell.text == h) {
                    headerRow = rk.first;
                    found = true;
                    break;
                }
            if (found)
                break;
        }
        if (found)
            break;
    }

    auto headerIt = rs.rows.find(headerRow);
    int maxHeaderCol = 0;
    if (headerIt != rs.rows.end() && !headerIt->second.cells.empty())
        maxHeaderCol = headerIt->second.cells.rbegin()->first;

    // 定位/新建目标列
    std::vector<int> targetCols(headers.size(), 0);
    for (std::size_t i = 0; i < headers.size(); ++i) {
        int col = -1;
        if (headerIt != rs.rows.end()) {
            for (const auto& ck : headerIt->second.cells)
                if (ck.second.value.isString &&
                    ck.second.value.text == headers[i]) {
                    col = ck.first;
                    break;
                }
        }
        if (col < 0) {
            col = ++maxHeaderCol;
            RawCell rc;
            rc.raw = "<c r=\"" + colName(col) + std::to_string(headerRow) +
                     "\" t=\"inlineStr\"><is><t>" + escapeXml(headers[i]) +
                     "</t></is></c>";
            rc.value.hasValue = true;
            rc.value.isString = true;
            rc.value.text = headers[i];
            rs.rows[headerRow].cells[col] = rc;
        }
        targetCols[i] = col;
    }

    // 公共起始行 = 各目标列已有数据末行最大值 + 1(下界为表头行+1)
    int startRow = headerRow + 1;
    for (int col : targetCols) {
        for (const auto& rk : rs.rows) {
            if (rk.first <= headerRow)
                continue;
            const auto it = rk.second.cells.find(col);
            if (it != rk.second.cells.end() && it->second.value.hasValue &&
                rk.first + 1 > startRow)
                startRow = rk.first + 1;
        }
    }

    // 追加写入数据
    const std::size_t rows = columns.empty() ? 0 : columns[0].size();
    for (std::size_t r = 0; r < rows; ++r) {
        const int rowNum = startRow + static_cast<int>(r);
        for (std::size_t i = 0; i < headers.size(); ++i) {
            if (std::isnan(columns[i][r]))
                continue;
            RawCell rc;
            rc.raw = "<c r=\"" + colName(targetCols[i]) +
                     std::to_string(rowNum) + "\"><v>" +
                     formatNumber(columns[i][r]) + "</v></c>";
            rc.value.hasValue = true;
            rc.value.isString = false;
            rc.value.number = columns[i][r];
            rs.rows[rowNum].cells[targetCols[i]] = rc;
        }
    }

    return rs.prefix + "<sheetData>" + serializeRows(rs.rows) +
           "</sheetData>" + rs.suffix;
}

// 在包的 workbook / rels / content-types 中登记一张新工作表
void registerSheet(std::vector<OutPart>& parts, const std::string& wb,
                   const std::string& rels, const std::string& ct,
                   const std::string& sheetName,
                   const std::string& sheetPath) {
    int maxSheetId = 0;
    {
        std::size_t q = 0;
        for (;;) {
            std::size_t s = wb.find("sheetId=\"", q);
            if (s == std::string::npos)
                break;
            s += 9;
            const std::size_t e = wb.find('"', s);
            try {
                maxSheetId =
                    std::max(maxSheetId, std::stoi(wb.substr(s, e - s)));
            } catch (...) {
            }
            q = e + 1;
        }
    }
    const int newSheetId = maxSheetId + 1;

    int maxRid = 0;
    for (const std::string& src : {wb, rels}) {
        std::size_t q = 0;
        for (;;) {
            const std::size_t r = src.find("rId", q);
            if (r == std::string::npos)
                break;
            std::size_t j = r + 3;
            int val = 0;
            bool any = false;
            while (j < src.size() &&
                   std::isdigit(static_cast<unsigned char>(src[j]))) {
                val = val * 10 + (src[j] - '0');
                ++j;
                any = true;
            }
            if (any) {
                maxRid = std::max(maxRid, val);
                q = j;
            } else {
                q = r + 3;
            }
        }
    }
    const std::string newRid = "rId" + std::to_string(maxRid + 1);

    std::string newWb = wb;
    const std::size_t sp = newWb.find("</sheets>");
    newWb.insert(sp, "<sheet name=\"" + escapeXml(sheetName) +
                         "\" sheetId=\"" + std::to_string(newSheetId) +
                         "\" r:id=\"" + newRid + "\"/>");
    replacePart(parts, "xl/workbook.xml",
                makeDeflated("xl/workbook.xml", newWb));

    std::string newRels = rels;
    const std::size_t rp = newRels.find("</Relationships>");
    newRels.insert(rp, "<Relationship Id=\"" + newRid +
                           "\" Type=\"http://schemas.openxmlformats.org/"
                           "officeDocument/2006/relationships/worksheet\" "
                           "Target=\"" +
                           sheetPath.substr(3) + "\"/>");
    replacePart(parts, "xl/_rels/workbook.xml.rels",
                makeDeflated("xl/_rels/workbook.xml.rels", newRels));

    if (!ct.empty()) {
        std::string newCt = ct;
        const std::size_t cp = newCt.find("</Types>");
        newCt.insert(cp, "<Override PartName=\"/" + sheetPath +
                             "\" ContentType=\"application/vnd."
                             "openxmlformats-officedocument.spreadsheetml."
                             "worksheet+xml\"/>");
        replacePart(parts, "[Content_Types].xml",
                    makeDeflated("[Content_Types].xml", newCt));
    }
}

}  // namespace

void readPointsFromXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::string& xHeader,
                        const std::string& yHeader,
                        std::vector<double>& x,
                        std::vector<double>& y) {
    const std::vector<uint8_t> data = readFile(path);
    const std::vector<ZipEntry> entries = readCentralDirectory(data);

    std::vector<std::string> shared;
    if (const ZipEntry* ss = findEntry(entries, "xl/sharedStrings.xml"))
        shared = parseSharedStrings(toString(extractEntry(data, *ss)));

    const std::string sheetPath = resolveSheetPath(data, entries, sheetName);
    const ZipEntry* sheet = findEntry(entries, sheetPath);
    if (sheet == nullptr)
        throw std::runtime_error("worksheet not found: " + sheetPath);

    const std::string xml = toString(extractEntry(data, *sheet));
    const std::map<int, std::map<int, Cell>> cells = parseSheet(xml, shared);

    int headerRow = -1;
    int xCol = -1;
    int yCol = -1;
    for (const auto& row : cells) {
        for (const auto& col : row.second) {
            const Cell& cell = col.second;
            if (!cell.isString)
                continue;
            if (cell.text == xHeader) {
                headerRow = row.first;
                xCol = col.first;
            } else if (cell.text == yHeader) {
                headerRow = row.first;
                yCol = col.first;
            }
        }
    }
    if (headerRow < 0 || xCol < 0 || yCol < 0)
        throw std::runtime_error("header not found: " + xHeader + " / " +
                                 yHeader);

    x.clear();
    y.clear();
    for (const auto& row : cells) {
        if (row.first <= headerRow)
            continue;
        const auto itx = row.second.find(xCol);
        const auto ity = row.second.find(yCol);
        if (itx == row.second.end() || ity == row.second.end())
            continue;
        if (itx->second.isString || ity->second.isString)
            continue;
        if (!itx->second.hasValue || !ity->second.hasValue)
            continue;
        x.push_back(itx->second.number);
        y.push_back(ity->second.number);
    }
}

void writeColumnsToXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::vector<std::string>& headers,
                        const std::vector<std::vector<double>>& columns,
                        bool overwrite, int startRow) {
    if (headers.size() != columns.size())
        throw std::runtime_error("headers and columns size mismatch");

    // 各数据列长度必须一致
    for (std::size_t i = 1; i < columns.size(); ++i)
        if (columns[i].size() != columns[0].size())
            throw std::runtime_error(
                "all data columns must have the same size");

    // 读取整个包：未修改的部件保留原始压缩字节，仅在需要时解压(W1)
    std::vector<OutPart> parts;
    std::vector<uint8_t> data;
    std::vector<ZipEntry> entries;
    std::string workbookText, relsText, ctText;
    std::string targetPath;
    bool readable = false;
    bool sheetExists = false;

    try {
        data = readFile(path);
        entries = readCentralDirectory(data);
        parts.reserve(entries.size() + 2);
        for (const auto& e : entries) {
            OutPart p;
            p.name = e.name;
            p.method = e.method;
            p.crc = e.crc;
            p.compSize = e.compSize;
            p.uncompSize = e.uncompSize;
            try {
                p.data = rawEntryBytes(data, e);  // 原样保留，不解压
            } catch (...) {
                // 本地头不可用时退回解压后按 stored 写
                p.data = extractEntry(data, e);
                p.method = 0;
                p.uncompSize = static_cast<uint32_t>(p.data.size());
                p.compSize = p.uncompSize;
                p.crc = crc32(p.data.data(), p.data.size());
            }
            if (e.name == "xl/workbook.xml")
                workbookText = partText(p);
            else if (e.name == "xl/_rels/workbook.xml.rels")
                relsText = partText(p);
            else if (e.name == "[Content_Types].xml")
                ctText = partText(p);
            parts.push_back(std::move(p));
        }
        readable = !workbookText.empty();
        if (readable) {
            try {
                targetPath = resolveSheetPath(data, entries, sheetName);
                sheetExists = true;
            } catch (...) {
                sheetExists = false;
            }
        }
    } catch (...) {
        parts.clear();
        readable = false;
    }

    if (!readable) {
        parts = buildEmptyPackage();
        for (const auto& p : parts) {
            if (p.name == "xl/workbook.xml")
                workbookText = partText(p);
            else if (p.name == "xl/_rels/workbook.xml.rels")
                relsText = partText(p);
            else if (p.name == "[Content_Types].xml")
                ctText = partText(p);
        }
        sheetExists = false;
    }

    // 共享字符串(解读已有表头/单元格需要)
    std::vector<std::string> shared;
    for (const auto& p : parts)
        if (p.name == "xl/sharedStrings.xml") {
            try {
                shared = parseSharedStrings(partText(p));
            } catch (...) {
            }
            break;
        }

    // 取目标工作表 XML(存在时)，在内存中增量修补(C1)
    std::string sheetXml;
    if (sheetExists) {
        try {
            for (const auto& p : parts)
                if (p.name == targetPath) {
                    sheetXml = partText(p);
                    break;
                }
        } catch (...) {
            sheetXml.clear();
        }
    }

    if (!sheetExists) {
        // 分配新工作表路径
        int maxSheet = 0;
        const std::string pre = "xl/worksheets/sheet";
        const std::string suf = ".xml";
        for (const auto& p : parts) {
            if (p.name.rfind(pre, 0) == 0 &&
                p.name.size() > pre.size() + suf.size() &&
                p.name.compare(p.name.size() - suf.size(), suf.size(),
                               suf) == 0) {
                const std::string num = p.name.substr(
                    pre.size(), p.name.size() - pre.size() - suf.size());
                bool digits = !num.empty();
                for (const char ch : num)
                    if (!(ch >= '0' && ch <= '9'))
                        digits = false;
                if (digits)
                    maxSheet = std::max(maxSheet, std::stoi(num));
            }
        }
        targetPath =
            "xl/worksheets/sheet" + std::to_string(maxSheet + 1) + ".xml";
    }

    const std::string patched =
        patchSheetXml(sheetXml, shared, headers, columns, overwrite, startRow);
    const OutPart sheetPart = makeDeflated(targetPath, patched);

    if (sheetExists) {
        replacePart(parts, targetPath, sheetPart);
    } else {
        parts.push_back(sheetPart);
        registerSheet(parts, workbookText, relsText, ctText, sheetName,
                      targetPath);
    }

    writeZip(path, parts);
}
