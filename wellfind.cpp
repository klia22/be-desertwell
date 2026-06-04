#include <bits/stdc++.h>
#include "well.h"
#include <atomic>
#include <chrono>
#include <csignal>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

using namespace std;
using namespace std::chrono;

static atomic<bool> stopRequested{false};
static atomic<uint64_t> foundCount{0};
static mutex outputMutex;

enum Corner : int {
    NW = 0,
    NE = 1,
    SW = 2,
    SE = 3,
    CORNER_COUNT = 4
};

static constexpr uint64_t TOTAL_UINT32 = 0x1'0000'0000ULL;
static constexpr uint32_t REPORT_INTERVAL_SECONDS = 5;
static constexpr int LOW12_BITS = 12;
static constexpr uint32_t LOW12_SIZE = 1u << LOW12_BITS;
static constexpr uint32_t LOW12_MASK = LOW12_SIZE - 1u;
static constexpr int LOW12_WORDS = LOW12_SIZE / 64;

using Low12Mask = array<uint64_t, LOW12_WORDS>;

static constexpr array<pair<int, int>, CORNER_COUNT> cornerOffsets = {
    pair{1, 1},   // NW
    pair{14, 1},  // NE
    pair{1, 14},  // SW
    pair{14, 14}  // SE
};

void handleSigint(int) {
    stopRequested.store(true, memory_order_relaxed);
}

inline void low12Set(Low12Mask &mask, uint32_t idx) noexcept {
    mask[idx >> 6] |= 1ull << (idx & 63);
}

inline uint32_t computeRegionSeedFromBase(uint32_t base) noexcept {
    return base ^ (DW::FEATURE_KEY + (base << 6) + (base >> 2) - 1640531527u);
}

inline bool low12Any(const Low12Mask &mask) noexcept {
    for (uint64_t w : mask) {
        if (w) return true;
    }
    return false;
}

inline Low12Mask low12And(const Low12Mask &a, const Low12Mask &b) noexcept {
    Low12Mask result{};
    for (int i = 0; i < LOW12_WORDS; ++i) {
        result[i] = a[i] & b[i];
    }
    return result;
}

inline int countTrailingZeros(uint64_t value) noexcept {
    return value ? __builtin_ctzll(value) : 64;
}

inline void printProgress(const string &phaseName, uint64_t processed, uint64_t total, double rate) {
    lock_guard<mutex> guard(outputMutex);
    double percent = total ? (100.0 * processed / double(total)) : 0.0;
    cerr << "[" << phaseName << "] "
         << fixed << setprecision(2) << percent << "% (" << processed << "/" << total << ")"
         << " rate=" << rate << "/s";
    if (phaseName == "search") {
        cerr << " found=" << foundCount.load(memory_order_relaxed);
    }
    cerr << '\n';
}

struct BestSolution {
    uint32_t seed;
    uint32_t xMul;
    uint32_t zMul;
    uint32_t baseO;
    int32_t chunkX;
    int32_t chunkZ;
    uint64_t distanceSquared;
};

inline uint64_t squaredDistance(int32_t x, int32_t z) noexcept {
    int64_t xx = int64_t(x);
    int64_t zz = int64_t(z);
    return uint64_t(xx * xx + zz * zz);
}

bool solveLinearDiophantine(int64_t a, int64_t b, int64_t c, int64_t &x0, int64_t &y0, int64_t &g) {
    auto extgcd = [&](auto self, int64_t aa, int64_t bb) -> pair<int64_t, int64_t> {
        if (bb == 0) return {1, 0};
        auto p = self(self, bb, aa % bb);
        return pair{p.second, p.first - (aa / bb) * p.second};
    };
    auto uv = extgcd(extgcd, a, b);
    g = a * uv.first + b * uv.second;
    if (g < 0) {
        g = -g;
        uv.first = -uv.first;
        uv.second = -uv.second;
    }
    if (c % g != 0) return false;
    x0 = uv.first * (c / g);
    y0 = uv.second * (c / g);
    return true;
}

BestSolution nearestSolution(uint32_t seed, uint32_t xMul, uint32_t zMul, uint32_t baseO) {
    int64_t cx0 = 0;
    int64_t cz0 = 0;
    int64_t g = 0;
    if (!solveLinearDiophantine(int64_t(xMul), int64_t(zMul), int64_t(int32_t(baseO)), cx0, cz0, g)) {
        return BestSolution{seed, xMul, zMul, baseO, 0, 0, UINT64_MAX};
    }
    int64_t stepX = int64_t(zMul) / g;
    int64_t stepZ = int64_t(xMul) / g;
    long double ideal = -static_cast<long double>(cx0) / static_cast<long double>(stepX);
    int64_t bestK = llround(ideal);
    int64_t bestCx = cx0 + bestK * stepX;
    int64_t bestCz = cz0 - bestK * stepZ;
    uint64_t bestDistance = squaredDistance(int32_t(bestCx), int32_t(bestCz));
    for (int delta = -2; delta <= 2; ++delta) {
        int64_t candidateK = bestK + delta;
        int64_t candidateCx = cx0 + candidateK * stepX;
        int64_t candidateCz = cz0 - candidateK * stepZ;
        uint64_t candidateDistance = squaredDistance(int32_t(candidateCx), int32_t(candidateCz));
        if (candidateDistance < bestDistance) {
            bestDistance = candidateDistance;
            bestCx = candidateCx;
            bestCz = candidateCz;
        }
    }
    return BestSolution{seed, xMul, zMul, baseO, int32_t(bestCx), int32_t(bestCz), bestDistance};
}

struct PhaseStatus {
    string name;
    uint64_t total;
    atomic<uint64_t> processed{};
};

void reporterThread(const PhaseStatus &status) {
    uint64_t lastProcessed = 0;
    auto lastTime = steady_clock::now();
    while (!stopRequested.load(memory_order_relaxed)) {
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
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    unsigned threadCount = thread::hardware_concurrency();
    if (threadCount == 0) threadCount = 1;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            threadCount = unsigned(stoi(argv[++i]));
        }
    }

    signal(SIGINT, handleSigint);

    vector<vector<uint32_t>> validBases(CORNER_COUNT);
    PhaseStatus baseStatus{"base-scan", TOTAL_UINT32};

    auto baseWorker = [&](unsigned tid) {
        uint64_t blockSize = (TOTAL_UINT32 + threadCount - 1) / threadCount;
        uint64_t start = uint64_t(tid) * blockSize;
        uint64_t end = min(start + blockSize, TOTAL_UINT32);
        array<vector<uint32_t>, CORNER_COUNT> localLists{};
        uint64_t localProcessed = 0;

        for (uint64_t value = start; value < end && !stopRequested.load(memory_order_relaxed); ++value) {
            uint32_t base = uint32_t(value);
            uint32_t regionSeed = computeRegionSeedFromBase(base);
            DW::RandWrapper rng(regionSeed);
            if (rng.nextInt<500>() != 0u) {
                ++localProcessed;
                if ((localProcessed & 0x3FFFFF) == 0) {
                    baseStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
                    localProcessed = 0;
                }
                continue;
            }
            uint32_t offZ = rng.nextInt<16>();
            uint32_t offX = rng.nextInt<16>();
            for (int corner = 0; corner < CORNER_COUNT; ++corner) {
                if (offX == uint32_t(cornerOffsets[corner].first) && offZ == uint32_t(cornerOffsets[corner].second)) {
                    localLists[corner].push_back(base);
                    break;
                }
            }
            ++localProcessed;
            if ((localProcessed & 0x3FFFFF) == 0) {
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
    for (unsigned t = 0; t < threadCount; ++t) {
        baseThreads.emplace_back(baseWorker, t);
    }
    thread baseReporter(reporterThread, cref(baseStatus));

    for (auto &t : baseThreads) t.join();
    stopRequested.store(false, memory_order_relaxed);
    if (baseReporter.joinable()) baseReporter.join();

    array<vector<uint32_t>, CORNER_COUNT> cornerSets;
    array<Low12Mask, CORNER_COUNT> low12Mask{};
    array<array<vector<uint32_t>, LOW12_SIZE>, CORNER_COUNT> low12Buckets;

    for (int corner = 0; corner < CORNER_COUNT; ++corner) {
        auto &list = validBases[corner];
        sort(list.begin(), list.end());
        list.erase(unique(list.begin(), list.end()), list.end());
        cornerSets[corner] = move(list);
        for (uint32_t value : cornerSets[corner]) {
            low12Set(low12Mask[corner], value & LOW12_MASK);
            low12Buckets[corner][value & LOW12_MASK].push_back(value);
        }
        cerr << "corner " << corner << " candidates=" << cornerSets[corner].size() << "\n";
        if (cornerSets[corner].empty()) {
            cerr << "ERROR: no candidates for corner " << corner << "\n";
            return 1;
        }
    }

    array<array<Low12Mask, LOW12_SIZE>, 3> xorLow12Masks{};
    array<array<uint8_t, LOW12_SIZE>, 3> low12Possible{};
    array<Low12Mask, 3> baseLow12Mask = {low12Mask[SW], low12Mask[NE], low12Mask[NW]};

    auto buildXorMasks = [&](int index, const Low12Mask &sourceMask) {
        for (uint32_t x = 0; x < LOW12_SIZE; ++x) {
            auto &mask = xorLow12Masks[index][x];
            for (uint32_t low12 = 0; low12 < LOW12_SIZE; ++low12) {
                uint32_t word = low12 >> 6;
                uint64_t bit = 1ull << (low12 & 63);
                if (sourceMask[word] & bit) {
                    uint32_t target = low12 ^ x;
                    mask[target >> 6] |= 1ull << (target & 63);
                }
            }
            low12Possible[index][x] = low12Any(low12And(low12Mask[SE], xorLow12Masks[index][x]));
        }
    };

    buildXorMasks(0, low12Mask[SW]); // xMul low12 from SW-SE
    buildXorMasks(1, low12Mask[NE]); // zMul low12 from NE-SE
    buildXorMasks(2, low12Mask[NW]); // xMul^zMul low12 from NW-SE

    // Recompute the possibility tables for the actual target sets.
    for (uint32_t x = 0; x < LOW12_SIZE; ++x) {
        low12Possible[0][x] = low12Any(low12And(low12Mask[SE], xorLow12Masks[0][x]));
        low12Possible[1][x] = low12Any(low12And(low12Mask[SE], xorLow12Masks[1][x]));
        low12Possible[2][x] = low12Any(low12And(low12Mask[SE], xorLow12Masks[2][x]));
    }

    PhaseStatus searchStatus{"search", TOTAL_UINT32};
    BestSolution best{0, 0, 0, 0, 0, 0, UINT64_MAX};
    mutex bestMutex;

    auto isValidBase = [&](uint32_t baseSE, uint32_t xMul, uint32_t zMul) {
        uint32_t baseSW = baseSE ^ xMul;
        uint32_t baseNE = baseSE ^ zMul;
        uint32_t baseNW = baseSE ^ xMul ^ zMul;
        if (!binary_search(cornerSets[SW].begin(), cornerSets[SW].end(), baseSW)) return false;
        if (!binary_search(cornerSets[NE].begin(), cornerSets[NE].end(), baseNE)) return false;
        if (!binary_search(cornerSets[NW].begin(), cornerSets[NW].end(), baseNW)) return false;
        return true;
    };

    auto searchWorker = [&](unsigned tid) {
        uint64_t blockSize = (TOTAL_UINT32 + threadCount - 1) / threadCount;
        uint64_t start = uint64_t(tid) * blockSize;
        uint64_t end = min(start + blockSize, TOTAL_UINT32);
        uint64_t localProcessed = 0;
        BestSolution localBest{0, 0, 0, 0, 0, 0, UINT64_MAX};

        for (uint64_t value = start; value < end && !stopRequested.load(memory_order_relaxed); ++value) {
            uint32_t seed = uint32_t(value);
            auto feat = DW::makeFeatureSeed(int64_t(seed));
            uint32_t xMul = feat.xMul;
            uint32_t zMul = feat.zMul;
            uint32_t xLo12 = xMul & LOW12_MASK;
            uint32_t zLo12 = zMul & LOW12_MASK;
            uint32_t xzLo12 = xLo12 ^ zLo12;

            if (!low12Possible[0][xLo12] || !low12Possible[1][zLo12] || !low12Possible[2][xzLo12]) {
                ++localProcessed;
                if ((localProcessed & 0x3FFFFF) == 0) {
                    searchStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
                    localProcessed = 0;
                }
                continue;
            }

            Low12Mask candidateMask{};
            for (int i = 0; i < LOW12_WORDS; ++i) {
                candidateMask[i] = low12Mask[SE][i] & xorLow12Masks[0][xLo12][i] & xorLow12Masks[1][zLo12][i] & xorLow12Masks[2][xzLo12][i];
            }
            if (!low12Any(candidateMask)) {
                ++localProcessed;
                if ((localProcessed & 0x3FFFFF) == 0) {
                    searchStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
                    localProcessed = 0;
                }
                continue;
            }

            for (int word = 0; word < LOW12_WORDS; ++word) {
                uint64_t w = candidateMask[word];
                while (w) {
                    int bit = countTrailingZeros(w);
                    w &= w - 1;
                    uint32_t low12 = (word << 6) + bit;
                    for (uint32_t baseSE : low12Buckets[SE][low12]) {
                        if (!isValidBase(baseSE, xMul, zMul)) continue;
                        uint32_t baseO = seed ^ baseSE;
                        BestSolution solution = nearestSolution(seed, xMul, zMul, baseO);
                        if (solution.distanceSquared < localBest.distanceSquared) {
                            localBest = solution;
                        }
                        foundCount.fetch_add(1, memory_order_relaxed);
                        break;
                    }
                    if (localBest.distanceSquared < UINT64_MAX) break;
                }
                if (localBest.distanceSquared < UINT64_MAX) break;
            }

            ++localProcessed;
            if ((localProcessed & 0x3FFFFF) == 0) {
                searchStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
                localProcessed = 0;
            }
        }
        searchStatus.processed.fetch_add(localProcessed, memory_order_relaxed);
        if (localBest.distanceSquared < UINT64_MAX) {
            lock_guard<mutex> guard(bestMutex);
            if (localBest.distanceSquared < best.distanceSquared) {
                best = localBest;
            }
        }
    };

    vector<thread> searchThreads;
    for (unsigned t = 0; t < threadCount; ++t) {
        searchThreads.emplace_back(searchWorker, t);
    }
    thread searchReporter(reporterThread, cref(searchStatus));

    for (auto &t : searchThreads) t.join();
    stopRequested.store(true, memory_order_relaxed);
    if (searchReporter.joinable()) searchReporter.join();

    if (best.distanceSquared < UINT64_MAX) {
        cerr << "FOUND seed=" << best.seed
             << " xMul=" << best.xMul
             << " zMul=" << best.zMul
             << " originChunk=(" << best.chunkX << "," << best.chunkZ << ")"
             << " distance^2=" << best.distanceSquared << "\n";
    } else {
        cerr << "No valid solution found in scanned range." << '\n';
    }
    return 0;
}
