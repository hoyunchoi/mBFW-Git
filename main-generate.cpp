#include <chrono>
#include <iostream>

#include "../library/stringFormat.hpp"

#include "generate.hpp"
#include "parameters.hpp"

int main(int argc, char* argv[]) {
    //* Get input parameters
    mBFW::generate::networkSize = std::stoi(argv[1]);
    mBFW::generate::acceptanceThreshold = std::stod(argv[2]);
    mBFW::generate::ensembleSize = std::stoul(argv[3]);
    mBFW::generate::coreNum = std::stoi(argv[4]);
    mBFW::generate::randomEngineSeed = -1;

    //* Check input network size and acceptance threshold
    if (mBFW::parameters::networkSizeList.find(mBFW::generate::networkSize) == mBFW::parameters::networkSizeList.end()) {
        std::cout << "WARNING: network size is not valid\n";
        return -1;
    }
    if (mBFW::parameters::acceptanceThresholdList.find(mBFW::generate::acceptanceThreshold) == mBFW::parameters::acceptanceThresholdList.end()) {
        std::cout << "WARNING: acceptance threshold is not valid\n";
        return -1;
    }

    //* Run mBFW
    auto start = std::chrono::system_clock::now();
    mBFW::generate::setParameters();
    mBFW::generate::run();
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    FILE* log = fopen("time.log", "a");
    fprintf(log, " %.6fs for N=%.1e, g=%.1f, ensemble=%d-%d\n", sec.count(), (double)mBFW::generate::networkSize, mBFW::generate::acceptanceThreshold, mBFW::generate::ensembleSize, mBFW::generate::coreNum);

    //* Save Data
    start = std::chrono::system_clock::now();
    mBFW::generate::save();
    sec = std::chrono::system_clock::now() - start;
    fprintf(log, " %0.6fs for saving\n", sec.count());

    return 0;
}
