#include <iostream>
#include <chrono>

#include "generate.hpp"

int main(int argc, char *argv[]){
    //* Get input parameters
    const int networkSize=std::stoul(argv[1]);
    const double acceptanceThreshold=std::stod(argv[2]);
    const int ensembleSize=std::stoul(argv[3]);
    const std::string machine=argv[4];
    const int coreNum=std::stoul(argv[5]);
    constexpr int randomEngineSeed = -1;

    //* Run mBFW
    auto start=std::chrono::system_clock::now();
    mBFW::generate::setParameters(networkSize, ensembleSize, acceptanceThreshold, coreNum, randomEngineSeed);
    mBFW::generate::run();
    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    FILE* log = fopen("log.txt", "a");
    fprintf(log, " %.6fs for N=%.1e, g=%.1f, ensemble=%d-%d at %s\n", sec.count(),(double)networkSize, acceptanceThreshold, ensembleSize, coreNum, machine.c_str());

    //* Save Data
    start = std::chrono::system_clock::now();
    mBFW::generate::save();
    sec=std::chrono::system_clock::now()-start;
    fprintf(log, " %0.6fs for saving\n", sec.count());

    return 0;
}
