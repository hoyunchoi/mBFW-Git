#include <iostream>
#include <chrono>

#include "data.hpp"

int main(int argc, char *argv[]){
    // const int networkSize=1000;
    // const double acceptanceThreshold = 0.5;
    const int networkSize=std::stoul(argv[1]);
    const double acceptanceThreshold=std::stod(argv[2]);
    const bool deletion = true;
    const double logBinDelta = 0.1;

    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, logBinDelta, deletion);
    mBFW::data::process();
    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}