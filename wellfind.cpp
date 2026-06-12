#include <bits/stdc++.h>
#include "well.h" //Used only for RNG
#include <atomic>
#include <chrono>
#include <csignal>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>
#include <numeric>
#include <immintrin.h>
#include <memory>

using namespace std;
using namespace std::chrono;

static atomic<bool> stopRequested{false};
static atomic<uint64_t> foundCount{0};
static atomic<uint64_t> candidateCount{0};
static mutex outputMutex;

enum Corner : int {
    NW = 0,
    NE = 1,
    SW = 2,
    SE = 3,
    CORNER_COUNT = 4
};

static constexpr uint64_t TOTAL_UINT32 = 0x1'0000'0000ULL;
static constexpr uint32_t REPORT_INTERVAL_SECONDS = 20;
static constexpr uint64_t FLUSH_INTERVAL = 20;

// Exact low22 filtering
static constexpr int LOWER_BITS = 22;
static constexpr uint32_t LOWER_SIZE = 1u << LOWER_BITS;
static constexpr uint32_t LOWER_MASK = LOWER_SIZE - 1u;

// MITM split for the low22 layer: 12 high bits + 10 low bits
static constexpr int MITM_LOW_BITS = 10;
static constexpr int MITM_HIGH_BITS = LOWER_BITS - MITM_LOW_BITS; // 12
static constexpr uint32_t MITM_LOW_SIZE = 1u << MITM_LOW_BITS;    // 1024
static constexpr uint32_t MITM_HIGH_SIZE = 1u << MITM_HIGH_BITS;  // 4096
static constexpr uint32_t MITM_LOW_MASK = MITM_LOW_SIZE - 1u;
static constexpr uint32_t MITM_HIGH_MASK = MITM_HIGH_SIZE - 1u;

void handleSigint(int) {
    stopRequested.store(true, memory_order_relaxed);
}

inline uint32_t computeRegionSeedFromBase(uint32_t base) noexcept {
    return base ^ (DW::FEATURE_KEY + (base << 6) + (base >> 2) - 1640531527u);
}

// Fast MT19937
static constexpr uint32_t MT_A = 1812433253u;
static constexpr uint32_t TWIST_B = 0x9908b0dfu;

inline uint32_t temper(uint32_t y) noexcept {
    y ^= (y >> 11);
    y ^= (y << 7) & 2636928640u;
    y ^= (y << 15) & 4022730752u;
    y ^= (y >> 18);
    return y;
}

inline uint32_t twistOnce(uint32_t a, uint32_t b, uint32_t c) noexcept {
    uint32_t y = (a & 0x80000000u) | (b & 0x7fffffffu);
    uint32_t m = c ^ (y >> 1);
    if (y & 1u) m ^= TWIST_B;
    return m;
}

inline DW::FeatureSeed fastMakeFeatureSeed(uint32_t seedLow32) noexcept {
    uint32_t mt0 = seedLow32;
    uint32_t mt1 = 0, mt2 = 0, mt397 = 0, mt398 = 0;

    uint32_t prev = seedLow32;
    for (uint32_t i = 1; i <= 398; ++i) {
        prev = MT_A * (prev ^ (prev >> 30)) + i;
        if (i == 1) mt1 = prev;
        else if (i == 2) mt2 = prev;
        else if (i == 397) mt397 = prev;
        else if (i == 398) mt398 = prev;
    }

    uint32_t raw0 = temper(twistOnce(mt0, mt1, mt397));
    uint32_t raw1 = temper(twistOnce(mt1, mt2, mt398));

    return DW::FeatureSeed{
        seedLow32,
        (raw0 >> 1) | 1u,
        (raw1 >> 1) | 1u,
        DW::FEATURE_KEY
    };
}

inline void printProgress(const string &phaseName, uint64_t processed, uint64_t total, double rate) {
    lock_guard<mutex> guard(outputMutex);
    double percent = total ? (100.0 * processed / double(total)) : 0.0;
    cerr << '[' << phaseName << "] "
         << fixed << setprecision(3) << percent << "% (" << processed << '/' << total << ')'
         << " rate=" << rate << "/s";
    if (phaseName == "search") {
        cerr << " found=" << foundCount.load(memory_order_relaxed)
             << " candidates=" << candidateCount.load(memory_order_relaxed);
    }
    cerr << '\n';
}

bool solveLinearDiophantine(int64_t a, int64_t b, int64_t c,
                            int64_t &x0, int64_t &y0, int64_t &g) {
    auto extgcd = [&](auto self, int64_t aa, int64_t bb) -> pair<int64_t, int64_t> {
        if (bb == 0) return {1, 0};
        auto p = self(self, bb, aa % bb);
        __int128 x = p.second;
        __int128 y = (__int128)p.first - (__int128)(aa / bb) * p.second;
        return { (int64_t)x, (int64_t)y };
    };

    g = std::gcd(a, b);
    if (g < 0) g = -g;
    if (g == 0) return c == 0;
    if (c % g != 0) return false;

    auto uv = extgcd(extgcd, a, b);
    __int128 scale = (__int128)c / g;
    x0 = (int64_t)((__int128)uv.first * scale);
    y0 = (int64_t)((__int128)uv.second * scale);
    return true;
}

struct BestSolution {
    uint32_t seed;
    uint32_t xMul;
    uint32_t zMul;
    uint32_t baseO;
    int32_t chunkX;
    int32_t chunkZ;
    uint64_t distance;
};

static inline __int128 iabs128(__int128 v) noexcept {
    return v < 0 ? -v : v;
}

static inline __int128 floorDiv128(__int128 a, __int128 b) noexcept {
    __int128 q = a / b;
    __int128 r = a % b;
    if (r != 0 && a < 0) --q;
    return q;
}

static inline __int128 ceilDiv128(__int128 a, __int128 b) noexcept {
    __int128 q = a / b;
    __int128 r = a % b;
    if (r != 0 && a > 0) ++q;
    return q;
}

BestSolution nearestSolution(uint32_t seed, uint32_t xMul, uint32_t zMul, uint32_t baseO) {
    int64_t cx0 = 0;
    int64_t cz0 = 0;
    int64_t g = 0;

    if (!solveLinearDiophantine(int64_t(xMul), int64_t(zMul), int64_t(int32_t(baseO)), cx0, cz0, g)) {
        return BestSolution{seed, xMul, zMul, baseO, 0, 0, UINT64_MAX};
    }
    if (g <= 0) {
        return BestSolution{seed, xMul, zMul, baseO, 0, 0, UINT64_MAX};
    }

    int64_t stepX = int64_t(zMul) / g;
    int64_t stepZ = int64_t(xMul) / g;

    __int128 ax0 = cx0;
    __int128 az0 = cz0;
    __int128 sx = stepX;
    __int128 sz = stepZ;

    auto feasible = [&](uint64_t D, __int128 &loK, __int128 &hiK) -> bool {
        __int128 d = (__int128)D;
        __int128 lo1 = ceilDiv128(-d - ax0, sx);
        __int128 hi1 = floorDiv128(d - ax0, sx);
        __int128 lo2 = ceilDiv128(az0 - d, sz);
        __int128 hi2 = floorDiv128(az0 + d, sz);
        loK = max(lo1, lo2);
        hiK = min(hi1, hi2);
        return loK <= hiK;
    };

    uint64_t upper = uint64_t(max(iabs128(ax0), iabs128(az0)));
    uint64_t lo = 0, hi = upper;

    while (lo < hi) {
        uint64_t mid = lo + ((hi - lo) >> 1);
        __int128 kLo = 0, kHi = 0;
        if (feasible(mid, kLo, kHi)) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    __int128 kLo = 0, kHi = 0;
    if (!feasible(lo, kLo, kHi)) {
        return BestSolution{seed, xMul, zMul, baseO, 0, 0, UINT64_MAX};
    }

    __int128 k = kLo;
    __int128 bestCx = ax0 + k * sx;
    __int128 bestCz = az0 - k * sz;

    return BestSolution{
        seed, xMul, zMul, baseO, int32_t(bestCx), int32_t(bestCz), lo
    };
}

struct PhaseStatus {
    string name;
    uint64_t total;
    atomic<uint64_t> processed{};
};

void reporterThread(const PhaseStatus &status) {
    uint64_t lastProcessed = 0;
    auto lastTime = steady_clock::now();

    while (!stopRequested.load(memory_order_relaxed) &&
           status.processed.load(memory_order_relaxed) < status.total) {
        this_thread::sleep_for(seconds(REPORT_INTERVAL_SECONDS));
        if (stopRequested.load(memory_order_relaxed)) break;

        auto now = steady_clock::now();
        double elapsed = duration<double>(now - lastTime).count();
        uint64_t current = status.processed.load(memory_order_relaxed);
        double rate = elapsed > 0.0 ? double(current - lastProcessed) / elapsed : 0.0;
        printProgress(status.name, current, status.total, rate);
        lastTime = now;
        lastProcessed = current;
    }

    uint64_t current = status.processed.load(memory_order_relaxed);
    printProgress(status.name, current, status.total, 0.0);
}

//SIMD optimization
struct alignas(64) RowMask1024 {
    uint64_t w[16]{};
};

struct ShiftParams {
    unsigned wordShift;
    unsigned bitShift;
    unsigned bitShift2;
};

static inline bool anyMask(const RowMask1024 &m) noexcept {
    const __m256i* v = (const __m256i*)m.w;
    __m256i o01 = _mm256_or_si256(v[0], v[1]);
    __m256i o23 = _mm256_or_si256(v[2], v[3]);
    __m256i final_or = _mm256_or_si256(o01, o23);
    return !_mm256_testz_si256(final_or, final_or);
}

static inline void andMask2(RowMask1024 &a, const RowMask1024 &b) noexcept {
    __m256i* va = (__m256i*)a.w;
    const __m256i* vb = (const __m256i*)b.w;
    va[0] = _mm256_and_si256(va[0], vb[0]);
    va[1] = _mm256_and_si256(va[1], vb[1]);
    va[2] = _mm256_and_si256(va[2], vb[2]);
    va[3] = _mm256_and_si256(va[3], vb[3]);
}

static inline void orMask2(RowMask1024 &a, const RowMask1024 &b) noexcept {
    __m256i* va = (__m256i*)a.w;
    const __m256i* vb = (const __m256i*)b.w;
    va[0] = _mm256_or_si256(va[0], vb[0]);
    va[1] = _mm256_or_si256(va[1], vb[1]);
    va[2] = _mm256_or_si256(va[2], vb[2]);
    va[3] = _mm256_or_si256(va[3], vb[3]);
}

static inline __m256i swapBits256(__m256i v, __m256i mask, int shift) noexcept {
    __m256i s = _mm256_srli_epi64(v, shift);
    __m256i t = _mm256_and_si256(_mm256_xor_si256(s, v), mask);
    return _mm256_xor_si256(_mm256_xor_si256(v, t), _mm256_slli_epi64(t, shift));
}

static inline void xorPermute1024AVX2(RowMask1024 &m, uint32_t x) noexcept {
    __m256i v0 = _mm256_loadu_si256((__m256i*)&m.w[0]);
    __m256i v1 = _mm256_loadu_si256((__m256i*)&m.w[4]);
    __m256i v2 = _mm256_loadu_si256((__m256i*)&m.w[8]);
    __m256i v3 = _mm256_loadu_si256((__m256i*)&m.w[12]);

    if (x & 1u) {
        __m256i mask = _mm256_set1_epi64x(0x5555555555555555ull);
        v0 = swapBits256(v0, mask, 1); v1 = swapBits256(v1, mask, 1);
        v2 = swapBits256(v2, mask, 1); v3 = swapBits256(v3, mask, 1);
    }
    if (x & 2u) {
        __m256i mask = _mm256_set1_epi64x(0x3333333333333333ull);
        v0 = swapBits256(v0, mask, 2); v1 = swapBits256(v1, mask, 2);
        v2 = swapBits256(v2, mask, 2); v3 = swapBits256(v3, mask, 2);
    }
    if (x & 4u) {
        __m256i mask = _mm256_set1_epi64x(0x0f0f0f0f0f0f0f0full);
        v0 = swapBits256(v0, mask, 4); v1 = swapBits256(v1, mask, 4);
        v2 = swapBits256(v2, mask, 4); v3 = swapBits256(v3, mask, 4);
    }
    if (x & 8u) {
        __m256i mask = _mm256_set1_epi64x(0x00ff00ff00ff00ffull);
        v0 = swapBits256(v0, mask, 8); v1 = swapBits256(v1, mask, 8);
        v2 = swapBits256(v2, mask, 8); v3 = swapBits256(v3, mask, 8);
    }
    if (x & 16u) {
        __m256i mask = _mm256_set1_epi64x(0x0000ffff0000ffffull);
        v0 = swapBits256(v0, mask, 16); v1 = swapBits256(v1, mask, 16);
        v2 = swapBits256(v2, mask, 16); v3 = swapBits256(v3, mask, 16);
    }
    if (x & 32u) {
        __m256i mask = _mm256_set1_epi64x(0x00000000ffffffffull);
        v0 = swapBits256(v0, mask, 32); v1 = swapBits256(v1, mask, 32);
        v2 = swapBits256(v2, mask, 32); v3 = swapBits256(v3, mask, 32);
    }

    if (x & 64u) {
        v0 = _mm256_permute4x64_epi64(v0, 0xB1); v1 = _mm256_permute4x64_epi64(v1, 0xB1);
        v2 = _mm256_permute4x64_epi64(v2, 0xB1); v3 = _mm256_permute4x64_epi64(v3, 0xB1);
    }
    if (x & 128u) {
        v0 = _mm256_permute4x64_epi64(v0, 0x4E); v1 = _mm256_permute4x64_epi64(v1, 0x4E);
        v2 = _mm256_permute4x64_epi64(v2, 0x4E); v3 = _mm256_permute4x64_epi64(v3, 0x4E);
    }
    if (x & 256u) {
        std::swap(v0, v1); std::swap(v2, v3);
    }
    if (x & 512u) {
        std::swap(v0, v2); std::swap(v1, v3);
    }

    _mm256_storeu_si256((__m256i*)&m.w[0], v0);
    _mm256_storeu_si256((__m256i*)&m.w[4], v1);
    _mm256_storeu_si256((__m256i*)&m.w[8], v2);
    _mm256_storeu_si256((__m256i*)&m.w[12], v3);
}

static inline void rotateRight1024_opt_inplace(RowMask1024 &m, const ShiftParams &sp) noexcept {
    if (sp.bitShift == 0u) {
        if (sp.wordShift == 0u) return;
        RowMask1024 out;
        #pragma GCC unroll 16
        for (unsigned i = 0; i < 16; ++i) {
            out.w[i] = m.w[(i + sp.wordShift) & 15u];
        }
        m = out;
        return;
    }
    RowMask1024 out;
    #pragma GCC unroll 16
    for (unsigned i = 0; i < 16; ++i) {
        uint64_t a = m.w[(i + sp.wordShift) & 15u];
        uint64_t b = m.w[(i + sp.wordShift + 1u) & 15u];
        out.w[i] = (a >> sp.bitShift) | (b << sp.bitShift2);
    }
    m = out;
}

static inline RowMask1024 makePrefixMask1024(uint32_t bits) noexcept {
    RowMask1024 m{};
    if (bits == 0) return m;
    uint32_t fullWords = bits / 64u;
    uint32_t rem = bits % 64u;
    for (uint32_t i = 0; i < fullWords; ++i) m.w[i] = ~0ull;
    if (rem != 0u && fullWords < 16u) m.w[fullWords] = (1ull << rem) - 1ull;
    return m;
}

static array<RowMask1024, MITM_LOW_SIZE> carryMask0;
static array<RowMask1024, MITM_LOW_SIZE> carryMask1;

static void initCarryMasks() {
    for (uint32_t d = 0; d < MITM_LOW_SIZE; ++d) {
        carryMask0[d] = makePrefixMask1024(MITM_LOW_SIZE - d);
        for (int i = 0; i < 16; ++i) {
            carryMask1[d].w[i] = ~carryMask0[d].w[i];
        }
    }
}

// Fast Lookup Bitsets and Data Arrays
static vector<uint64_t> swBitset;
static vector<uint64_t> neBitset;
static vector<uint64_t> nwBitset;
static vector<uint32_t> seBucketStart;
static vector<uint16_t> seBucketCount;
static vector<uint32_t> seValuesFlat;

static array<vector<uint32_t>, CORNER_COUNT> cornerSets;
static array<vector<RowMask1024>, CORNER_COUNT> cornerRows;

static void buildCornerArtifacts(int corner) {
    auto &set = cornerSets[corner];
    auto &rows = cornerRows[corner];
    rows.assign(MITM_HIGH_SIZE, RowMask1024{});

    vector<uint32_t> sortedByLow22 = set;
    sort(sortedByLow22.begin(), sortedByLow22.end(), [](uint32_t a, uint32_t b) {
        uint32_t lowA = a & LOWER_MASK;
        uint32_t lowB = b & LOWER_MASK;
        if (lowA != lowB) return lowA < lowB;
        return a < b;
    });

    uint32_t prevLow22 = UINT32_MAX;
    for (uint32_t v : sortedByLow22) {
        uint32_t low22 = v & LOWER_MASK;
        if (low22 != prevLow22) {
            uint32_t row = low22 >> MITM_LOW_BITS;
            uint32_t bit = low22 & MITM_LOW_MASK;
            rows[row].w[bit >> 6] |= 1ull << (bit & 63u);
            prevLow22 = low22;
        }
    }

    if (corner == SW) {
        swBitset.assign(TOTAL_UINT32 / 64, 0);
        for (uint32_t v : set) swBitset[v >> 6] |= (1ULL << (v & 63u));
    } else if (corner == NE) {
        neBitset.assign(TOTAL_UINT32 / 64, 0);
        for (uint32_t v : set) neBitset[v >> 6] |= (1ULL << (v & 63u));
    } else if (corner == NW) {
        nwBitset.assign(TOTAL_UINT32 / 64, 0);
        for (uint32_t v : set) nwBitset[v >> 6] |= (1ULL << (v & 63u));
    } else if (corner == SE) {
        seBucketStart.assign(LOWER_SIZE, 0);
        seBucketCount.assign(LOWER_SIZE, 0);
        seValuesFlat = sortedByLow22; 
        
        for (size_t i = 0; i < seValuesFlat.size(); ) {
            uint32_t low22 = seValuesFlat[i] & LOWER_MASK;
            size_t j = i;
            while (j < seValuesFlat.size() && (seValuesFlat[j] & LOWER_MASK) == low22) ++j;
            seBucketStart[low22] = uint32_t(i);
            seBucketCount[low22] = uint16_t(j - i);
            i = j;
        }
    }
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    unsigned threadCount = thread::hardware_concurrency();
    if (threadCount == 0) threadCount = 1;
    uint64_t limitTotal = 1ull << 32;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            threadCount = unsigned(stoi(argv[++i]));
        } else if (arg == "--limit" && i + 1 < argc) {
            limitTotal = stoull(argv[++i]);
            if (limitTotal == 0) limitTotal = TOTAL_UINT32;
        }
    }

    signal(SIGINT, handleSigint);
    initCarryMasks();

    vector<vector<uint32_t>> validBases(CORNER_COUNT);
    PhaseStatus baseStatus{"base-scan", limitTotal};

    auto baseWorker = [&](unsigned tid) {
        uint64_t blockSize = (limitTotal + threadCount - 1) / threadCount;
        uint64_t start = uint64_t(tid) * blockSize;
        uint64_t end = min(start + blockSize, limitTotal);

        array<vector<uint32_t>, CORNER_COUNT> localLists{};
        uint64_t localProcessed = 0;

        for (uint64_t value = start; value < end && !stopRequested.load(memory_order_relaxed); ++value) {
            uint32_t base = uint32_t(value);
            uint32_t regionSeed = computeRegionSeedFromBase(base);
            DW::RandWrapper rng(regionSeed);

            if (rng.nextInt<500>() != 0u) {
                ++localProcessed;
                if ((localProcessed & (FLUSH_INTERVAL - 1)) == 0) {
                    baseStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
                    localProcessed = 0;
                }
                continue;
            }

            uint32_t offZ = rng.nextInt<16>();
            uint32_t offX = rng.nextInt<16>();

            auto inLowRange = [](uint32_t v) noexcept { return v <= 4u; };
            auto inHighRange = [](uint32_t v) noexcept { return v >= 11u && v <= 15u; };

            bool xLow = inLowRange(offX);
            bool xHigh = inHighRange(offX);
            bool zLow = inLowRange(offZ);
            bool zHigh = inHighRange(offZ);

            for (int corner = 0; corner < CORNER_COUNT; ++corner) {
                bool match = false;
                switch (corner) {
                    case NW: match = xLow && zLow; break;
                    case NE: match = xHigh && zLow; break;
                    case SW: match = xLow && zHigh; break;
                    case SE: match = xHigh && zHigh; break;
                }
                if (match) {
                    localLists[corner].push_back(base);
                    break;
                }
            }

            ++localProcessed;
            if ((localProcessed & (FLUSH_INTERVAL - 1)) == 0) {
                baseStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
                localProcessed = 0;
            }
        }

        baseStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
        lock_guard<mutex> guard(outputMutex);
        for (int corner = 0; corner < CORNER_COUNT; ++corner) {
            auto &globalList = validBases[corner];
            auto &localList = localLists[corner];
            globalList.insert(globalList.end(), localList.begin(), localList.end());
        }
    };

    vector<thread> baseThreads;
    baseThreads.reserve(threadCount);
    for (unsigned t = 0; t < threadCount; ++t) {
        baseThreads.emplace_back(baseWorker, t);
    }
    thread baseReporter(reporterThread, cref(baseStatus));

    for (auto &t : baseThreads) t.join();
    stopRequested.store(false, memory_order_relaxed);
    if (baseReporter.joinable()) baseReporter.join();

    cerr << "base phase complete\n";

    for (int corner = 0; corner < CORNER_COUNT; ++corner) {
        auto &list = validBases[corner];
        sort(list.begin(), list.end());
        list.erase(unique(list.begin(), list.end()), list.end());
        cornerSets[corner] = move(list);

        buildCornerArtifacts(corner);

        cerr << "corner " << corner << " candidates=" << cornerSets[corner].size() << "\n";
        if (cornerSets[corner].empty()) {
            cerr << "ERROR: no candidates for corner " << corner << "\n";
            return 1;
        }
    }

    PhaseStatus searchStatus{"search", limitTotal};
    BestSolution best{0, 0, 0, 0, 0, 0, UINT64_MAX};
    mutex bestMutex;

    auto searchWorker = [&](unsigned tid) {
        uint64_t blockSize = (limitTotal + threadCount - 1) / threadCount;
        blockSize = (blockSize + 1023u) & ~1023ULL; 
        uint64_t start = min(uint64_t(tid) * blockSize, limitTotal);
        uint64_t end = min(start + blockSize, limitTotal);

        uint64_t localProcessed = 0;
        BestSolution localBest{0, 0, 0, 0, 0, 0, UINT64_MAX};

        // Preallocate masks securely on heap per thread to avoid stack overflows
        auto pre_seRows = make_unique<RowMask1024[]>(MITM_HIGH_SIZE);
        auto pre_swRows = make_unique<RowMask1024[]>(MITM_HIGH_SIZE);
        auto pre_neRows = make_unique<RowMask1024[]>(MITM_HIGH_SIZE);
        auto pre_nwRows = make_unique<RowMask1024[]>(MITM_HIGH_SIZE);

        const auto *seRows = cornerRows[SE].data();
        const auto *swRows = cornerRows[SW].data();
        const auto *neRows = cornerRows[NE].data();
        const auto *nwRows = cornerRows[NW].data();

        // Alter search order so seedLow10 (inner cycle mapping) acts as outer block
        for (uint32_t sl10 = 0; sl10 < 1024; ++sl10) {
            if (stopRequested.load(memory_order_relaxed)) break;

            for (uint32_t i = 0; i < MITM_HIGH_SIZE; ++i) {
                pre_seRows[i] = seRows[i]; xorPermute1024AVX2(pre_seRows[i], sl10);
                pre_swRows[i] = swRows[i]; xorPermute1024AVX2(pre_swRows[i], sl10);
                pre_neRows[i] = neRows[i]; xorPermute1024AVX2(pre_neRows[i], sl10);
                pre_nwRows[i] = nwRows[i]; xorPermute1024AVX2(pre_nwRows[i], sl10);
            }

            for (uint64_t base = start; base < end; base += 1024) {
                if (stopRequested.load(memory_order_relaxed)) break;
                uint64_t currentSeedVal = base + sl10;
                if (currentSeedVal >= end) continue;

                uint32_t seed = uint32_t(currentSeedVal);
                auto feat = fastMakeFeatureSeed(seed);

                uint32_t xMul = feat.xMul;
                uint32_t zMul = feat.zMul;

                uint32_t seedLow24 = seed & LOWER_MASK;
                uint32_t seedHigh14 = (seedLow24 >> MITM_LOW_BITS) & MITM_HIGH_MASK;

                uint32_t xLow24 = xMul & LOWER_MASK;
                uint32_t zLow24 = zMul & LOWER_MASK;
                uint32_t xzLow24 = (xLow24 + zLow24) & LOWER_MASK;

                uint32_t xLow10 = xLow24 & MITM_LOW_MASK;
                uint32_t zLow10 = zLow24 & MITM_LOW_MASK;
                uint32_t xzLow10 = xzLow24 & MITM_LOW_MASK;

                uint32_t xHigh14 = (xLow24 >> MITM_LOW_BITS) & MITM_HIGH_MASK;
                uint32_t zHigh14 = (zLow24 >> MITM_LOW_BITS) & MITM_HIGH_MASK;
                uint32_t xzHigh14 = (xzLow24 >> MITM_LOW_BITS) & MITM_HIGH_MASK;

                auto getShiftParams = [](unsigned r) -> ShiftParams {
                    r &= 1023u;
                    return { r >> 6, r & 63u, 64u - (r & 63u) };
                };
                ShiftParams spSW = getShiftParams(xLow10);
                ShiftParams spNE = getShiftParams(zLow10);
                ShiftParams spNW = getShiftParams(xzLow10);

                for (uint32_t row0 = 0; row0 < MITM_HIGH_SIZE; ++row0) {
                    uint32_t uHigh14 = seedHigh14 ^ row0;
                    RowMask1024 mask = pre_seRows[row0];

                    // Process SW exactly equivalent to before but in-lined to easily abort sequentially
                    uint32_t sumSW = (uHigh14 + xHigh14) & MITM_HIGH_MASK;
                    uint32_t rowSW = seedHigh14 ^ sumSW;
                    RowMask1024 extraSW = pre_swRows[rowSW];
                    rotateRight1024_opt_inplace(extraSW, spSW);
                    andMask2(extraSW, carryMask0[xLow10]);

                    if (xLow10 != 0u) {
                        uint32_t rowSW1 = seedHigh14 ^ ((sumSW + 1u) & MITM_HIGH_MASK);
                        RowMask1024 extraSW1 = pre_swRows[rowSW1];
                        rotateRight1024_opt_inplace(extraSW1, spSW);
                        andMask2(extraSW1, carryMask1[xLow10]);
                        orMask2(extraSW, extraSW1);
                    }
                    andMask2(mask, extraSW);
                    if (!anyMask(mask)) continue;

                    // NE Phase
                    uint32_t sumNE = (uHigh14 + zHigh14) & MITM_HIGH_MASK;
                    uint32_t rowNE = seedHigh14 ^ sumNE;
                    RowMask1024 extraNE = pre_neRows[rowNE];
                    rotateRight1024_opt_inplace(extraNE, spNE);
                    andMask2(extraNE, carryMask0[zLow10]);

                    if (zLow10 != 0u) {
                        uint32_t rowNE1 = seedHigh14 ^ ((sumNE + 1u) & MITM_HIGH_MASK);
                        RowMask1024 extraNE1 = pre_neRows[rowNE1];
                        rotateRight1024_opt_inplace(extraNE1, spNE);
                        andMask2(extraNE1, carryMask1[zLow10]);
                        orMask2(extraNE, extraNE1);
                    }
                    andMask2(mask, extraNE);
                    if (!anyMask(mask)) continue;

                    // NW Phase
                    uint32_t sumNW = (uHigh14 + xzHigh14) & MITM_HIGH_MASK;
                    uint32_t rowNW = seedHigh14 ^ sumNW;
                    RowMask1024 extraNW = pre_nwRows[rowNW];
                    rotateRight1024_opt_inplace(extraNW, spNW);
                    andMask2(extraNW, carryMask0[xzLow10]);

                    if (xzLow10 != 0u) {
                        uint32_t rowNW1 = seedHigh14 ^ ((sumNW + 1u) & MITM_HIGH_MASK);
                        RowMask1024 extraNW1 = pre_nwRows[rowNW1];
                        rotateRight1024_opt_inplace(extraNW1, spNW);
                        andMask2(extraNW1, carryMask1[xzLow10]);
                        orMask2(extraNW, extraNW1);
                    }
                    andMask2(mask, extraNW);
                    if (!anyMask(mask)) continue;

                    for (int word = 0; word < 16; ++word) {
                        uint64_t bits = mask.w[word];
                        while (bits) {
                            unsigned bit = unsigned(__builtin_ctzll(bits));
                            bits &= (bits - 1);

                            uint32_t uLow10 = uint32_t(word * 64 + bit);
                            uint32_t u24 = (uHigh14 << MITM_LOW_BITS) | uLow10;
                            uint32_t baseSELow24 = seedLow24 ^ u24;
                            
                            uint32_t startIdx = seBucketStart[baseSELow24];
                            uint32_t count = seBucketCount[baseSELow24];

                            for (uint32_t i = 0; i < count; ++i) {
                                if ((i & 1023u) == 0u && stopRequested.load(memory_order_relaxed)) break;

                                uint32_t baseSE = seValuesFlat[startIdx + i];
                                uint32_t cSE = seed ^ baseSE;

                                uint32_t baseSW = seed ^ (cSE + xMul);
                                if (!((swBitset[baseSW >> 6] >> (baseSW & 63)) & 1)) continue;

                                uint32_t baseNE = seed ^ (cSE + zMul);
                                if (!((neBitset[baseNE >> 6] >> (baseNE & 63)) & 1)) continue;

                                uint32_t baseNW = seed ^ (cSE + xMul + zMul);
                                if (!((nwBitset[baseNW >> 6] >> (baseNW & 63)) & 1)) continue;

                                candidateCount.fetch_add(1, memory_order_relaxed);
                                foundCount.fetch_add(1, memory_order_relaxed);

                                uint32_t baseO = cSE;
                                BestSolution solution = nearestSolution(seed, xMul, zMul, baseO);
                                if (solution.distance < localBest.distance) {
                                    localBest = solution;
                                }
                            }
                        }
                    }
                }

                ++localProcessed;
                if ((localProcessed & (FLUSH_INTERVAL - 1)) == 0) {
                    searchStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
                    localProcessed = 0;
                }
            }
        }
        
        searchStatus.processed.fetch_add(localProcessed, memory_order_relaxed);

        if (localBest.distance < UINT64_MAX) {
            lock_guard<mutex> guard(bestMutex);
            if (localBest.distance < best.distance) {
                best = localBest;
            }
        }
    };

    vector<thread> searchThreads;
    searchThreads.reserve(threadCount);
    for (unsigned t = 0; t < threadCount; ++t) {
        searchThreads.emplace_back(searchWorker, t);
    }
    thread searchReporter(reporterThread, cref(searchStatus));

    for (auto &t : searchThreads) t.join();
    stopRequested.store(true, memory_order_relaxed);
    if (searchReporter.joinable()) searchReporter.join();

    if (best.distance < UINT64_MAX) {
        cerr << "FOUND seed=" << best.seed
             << " xMul=" << best.xMul
             << " zMul=" << best.zMul
             << " originChunk=(" << best.chunkX << ',' << best.chunkZ << ')'
             << " distance=" << best.distance << '\n';
    } else {
        cerr << "No valid solution found in scanned range.\n";
    }

    int a;
    cin >> a;
    return 0;
}
