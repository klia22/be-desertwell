#include "well.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

static std::atomic<bool> stopRequested{false};
static std::mutex outputMutex;

void handle_sigint(int) {
    stopRequested.store(true, std::memory_order_relaxed);
}

inline bool has2x2(std::uint32_t top, std::uint32_t bottom) noexcept {
    std::uint32_t common = top & bottom;
    return (common & (common >> 1u)) != 0u;
}

inline std::uint32_t firstMTRandomInt(std::uint32_t seed) noexcept {
    std::uint32_t mt0 = seed;
    std::uint32_t mt1 = std::uint32_t(DW::imul32(1812433253u, mt0 ^ (mt0 >> 30u)) + 1u);
    std::uint32_t prev = mt1;
    std::uint32_t mt397 = 0;
    for (int i = 2; i <= 397; ++i) {
        std::uint32_t curr = std::uint32_t(DW::imul32(1812433253u, prev ^ (prev >> 30u)) + std::uint32_t(i));
        if (i == 397) mt397 = curr;
        prev = curr;
    }
    std::uint32_t y = (mt0 & 0x80000000u) | (mt1 & 0x7fffffffu);
    std::uint32_t x = mt397 ^ (y >> 1u);
    if (y & 1u) x ^= 0x9908b0dfu;
    x ^= (x >> 11u);
    x ^= (x << 7u) & 2636928640u;
    x ^= (x << 15u) & 4022730752u;
    x ^= (x >> 18u);
    return x;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;

    uint64_t offset = 0;
    uint64_t totalSeeds = 0x1'0000'0000ULL;
    bool printSeeds = true;
    bool progress = true;

    constexpr int W = 10;
    constexpr int H = 10;
    constexpr int minX = -5;
    constexpr int minZ = -5;

    std::atomic<uint64_t> processed{0};
    std::atomic<uint64_t> found{0};

    std::signal(SIGINT, handle_sigint);

    auto worker = [&](unsigned tid) {
        uint64_t localProcessed = 0;
        uint64_t localFound = 0;
        std::vector<std::uint32_t> localMatches;

        const uint64_t seedsPerThread = (totalSeeds + threads - 1) / threads;
        const uint64_t start = offset + uint64_t(tid) * seedsPerThread;
        const uint64_t end = std::min(start + seedsPerThread, totalSeeds);

        for (uint64_t s = start; s < end; ++s) {
            if (stopRequested.load(std::memory_order_relaxed)) break;
            std::uint32_t seed = std::uint32_t(s);

            auto feat = DW::makeFeatureSeed(std::int64_t(seed));
            std::uint32_t prevRow = 0;
            bool matched = false;

            const std::uint32_t seedLow = feat.seedLow32;
            const std::uint32_t xMul = feat.xMul;
            const std::uint32_t zMul = feat.zMul;
            const std::uint32_t fKey = feat.featureKey;

            for (int z = 0; z < H; ++z) {
                int cz = minZ + z;
                std::uint32_t row = 0;

                std::uint32_t combined = DW::imul32(std::uint32_t(minX), xMul)
                                        + DW::imul32(std::uint32_t(cz), zMul);

                for (int x = 0; x < W; ++x) {
                    std::uint32_t base = seedLow ^ combined;
                    std::uint32_t region = std::uint32_t(base ^ (fKey + (base << 6) + (base >> 2) - 1640531527u));
                    bool w = (firstMTRandomInt(region) % 500u) == 0u;
                    row |= (std::uint32_t)w << x;
                    combined += xMul;
                }

                if (z != 0 && has2x2(prevRow, row)) {
                    matched = true;
                    break;
                }
                prevRow = row;
            }

            if (matched) {
                ++localFound;
                found.fetch_add(1, std::memory_order_relaxed);
                if (printSeeds) localMatches.push_back(seed);
            }

            ++localProcessed;
            if ((localProcessed & 0x3FF) == 0) {
                processed.fetch_add(localProcessed, std::memory_order_relaxed);
                localProcessed = 0;
            }
        }

        if (localProcessed) processed.fetch_add(localProcessed, std::memory_order_relaxed);

        if (!localMatches.empty()) {
            std::lock_guard<std::mutex> guard(outputMutex);
            for (auto v : localMatches) std::cout << v << '\n';
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned i = 0; i < threads; ++i) {
        pool.emplace_back(worker, i);
    }

    std::thread reporterThread;
    if (progress) {
        reporterThread = std::thread([&]() {
            using clock = std::chrono::steady_clock;
            auto lastTime = clock::now();
            uint64_t lastProcessed = 0;
            while (!stopRequested.load(std::memory_order_relaxed)) {
                std::this_thread::sleep_for(std::chrono::seconds(5));
                auto now = clock::now();
                uint64_t currentProcessed = processed.load(std::memory_order_relaxed);
                double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now - lastTime).count();
                uint64_t rate = elapsed > 0.0 ? uint64_t((currentProcessed - lastProcessed) / elapsed) : 0;
                lastTime = now;
                lastProcessed = currentProcessed;
                std::cerr << "processed=" << currentProcessed << " rate=" << rate << "/s found=" << found.load() << " threads=" << threads << "\n";
            }
        });
    }

    for (auto &t : pool) t.join();
    stopRequested.store(true);
    if (reporterThread.joinable()) reporterThread.join();
    std::cerr << "done processed=" << processed.load() << " found=" << found.load() << "\n";
    int a;
    std::cin>>a;
    return 0;
}
