#include <chrono>
#include <iostream>

#include "generate.hpp"
#include "parameters.hpp"

int main(int argc, char* argv[]) {
    //* Get input parameters
    const int networkSize = std::stoul(argv[1]);
    const double acceptanceThreshold = std::stod(argv[2]);
    const int ensembleSize = std::stoul(argv[3]);
    const int coreNum = std::stoul(argv[4]);
    constexpr int randomEngineSeed = 0;

    //* Check input network size and acceptance threshold
    if (mBFW::parameters::networkSizeList.find(networkSize) == mBFW::parameters::networkSizeList.end()) {
        std::cout << "WARNING: network size is not valid\n";
        return -1;
    }
    if (mBFW::parameters::acceptanceThresholdList.find(acceptanceThreshold) == mBFW::parameters::acceptanceThresholdList.end()) {
        std::cout << "WARNING: acceptance threshold is not valid\n";
        return -1;
    }

    //* Run mBFW
    auto start = std::chrono::system_clock::now();
    mBFW::generate::setParameters(networkSize, ensembleSize, acceptanceThreshold, coreNum, randomEngineSeed);
    mBFW::generate::run();
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    FILE* log = fopen("log.txt", "a");
    fprintf(log, " %.6fs for N=%.1e, g=%.1f, ensemble=%d-%d\n", sec.count(), (double)networkSize, acceptanceThreshold, ensembleSize, coreNum);

    //* Save Data
    start = std::chrono::system_clock::now();
    mBFW::generate::save();
    sec = std::chrono::system_clock::now() - start;
    fprintf(log, " %0.6fs for saving\n", sec.count());

    return 0;
}
