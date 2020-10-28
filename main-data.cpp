#include <iostream>
#include <chrono>

#include "data.hpp"

int main(){
    const int networkSize=1000;
    const double acceptanceThreshold = 0.5;
    const bool deletion = true;

    const double logBinDelta = 0.1;
    std::vector<int> ensembleList(20,10000);
    // ensembleList[0] = 5000;

    //* Determine which observables to calculate

    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, deletion);
    mBFW::data::process();
    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}