#include <bits/stdc++.h>
#include "well.h"
using namespace std;
using namespace DW;
int main() {
    uint32_t worldSeed = 2589;
    uint32_t chunkX = 273983839;
    uint32_t chunkZ = -300695624;
    for(int ax = chunkX - 1; ax <= chunkX + 1; ++ax) {
        for(int az = chunkZ - 1; az <= chunkZ + 1; ++az) {
            auto result = findwell(worldSeed, ax, az);
            if(result.hasWell) {
                cout << "Found well at world coordinates (" << result.worldX << ", " << result.worldZ << ") "<< "Chunk (" << ax << ", " << az << ")\n";
            } else {
                cout << "No well at chunk (" << ax << ", " << az << ")\n";
            }
        }
    }
}