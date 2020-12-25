// #include <iostream>
#include <chrono>

#include "data.hpp"

int main(int argc, char *argv[]){
    const double acceptanceThreshold=std::stod(argv[1]);
    const int networkSize=std::stoul(argv[2]);
    const bool deletion = true;
    const double logBinDelta = 0.1;

    //* Check list of each observables
    std::map<std::string, bool> checkList;
    checkList["orderParameter"] = false;
    checkList["orderParameter_trial"] = false;
    checkList["meanClusterSize"] = false;
    checkList["meanClusterSize_trial"] = false;
    checkList["interEventTime"] = false;
    checkList["orderParameterVariance"] = false;
    checkList["orderParameterVariance_trial"] = false;
    checkList["clusterSizeDist"] = false;
    checkList["clusterSizeDist_exact"] = false;
    checkList["clusterSizeDist_time"] = false;
    checkList["interEventTimeDist_time"] = true;
    checkList["ageDist_time"] = false;
    checkList["deltaUpperBoundDist_time"] = false;
    checkList["orderParameterDist"] = false;


    //* Run data process
    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, logBinDelta, deletion);
    mBFW::data::process(checkList);
    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}