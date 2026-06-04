#include "well.h"
#include <iostream>
#include <vector>
#include <tuple>

int main() {
    for(int i = 0; i < 1000000; i++) {
        if(DW::isWell(i,0,0) && DW::isWell(i,1,0)) std::cout << i << " has wells at (0,0) and (1,0)" << std::endl;
    }
}