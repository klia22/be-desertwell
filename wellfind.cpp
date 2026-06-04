#include <bits/stdc++.h>
#include "well.h"
#include <mutex>
#include <thread>
#include <vector>
using namespace std;

mutex cout_mutex;
vector<uint32_t> validregions[2][2]{};
void worker1(uint32_t start, uint32_t end) {
    for (uint32_t i = start; i < end; i++) {
        DW::RandWrapper rng(i);
        if (i % 100000000 == 0){
            {
                lock_guard<mutex> lock(cout_mutex);
                cout << i << " seeds checked\n";
            }
        }
        if (rng.nextInt<500>() != 0u) {
            continue;
        }
        int oz = int(rng.nextInt<16>());
        int ox = int(rng.nextInt<16>());
        if ((1 <= oz && oz <= 14) || (1 <= ox && ox <= 14)) {
            continue;
        }
        {
            lock_guard<mutex> lock(cout_mutex);
            validregions[ox/15][oz/15].push_back(i);
        }
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    uint32_t num_threads = thread::hardware_concurrency();
    uint32_t total_seeds = 0x1'0000'0000ULL - 1;
    uint32_t chunk_size = total_seeds / num_threads;
    
    vector<thread> threads;
    for (uint32_t t = 0; t < num_threads; t++) {
        uint32_t start = t * chunk_size;
        uint32_t end = (t == num_threads - 1) ? total_seeds : (t + 1) * chunk_size;
        threads.emplace_back(worker1, start, end);
    }
    for (auto& t : threads) {
        t.join();
    }
    return 0;
}