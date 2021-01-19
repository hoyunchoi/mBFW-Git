// #include <iostream>
#include <chrono>

#include "parameters.hpp"
#include "data.hpp"

int main(int argc, char *argv[]){
    const double acceptanceThreshold=std::stod(argv[1]);
    const int networkSize=std::stoul(argv[2]);
    const bool deletion = true;
    const double logBinDelta = 0.1;

    //* Check input network size and acceptance threshold
    if (mBFW::parameters::networkSizeList.find(networkSize) == mBFW::parameters::networkSizeList.end()){
        std::cout << "WARNING: network size is not valid\n";
        return -1;
    }
    if (mBFW::parameters::acceptanceThresholdList.find(acceptanceThreshold) == mBFW::parameters::acceptanceThresholdList.end()){
        std::cout << "WARNING: acceptance threshold is not valid\n";
        return -1;
    }

    //* Check list of each observables
    std::map<std::string, bool> checkList;
    checkList["orderParameter"] = true;
    checkList["orderParameter_trial"] = true;
    checkList["meanClusterSize"] = true;
    checkList["meanClusterSize_trial"] = true;
    checkList["interEventTime"] = true;
    checkList["orderParameterVariance"] = true;
    checkList["orderParameterVariance_trial"] = true;
    checkList["clusterSizeDist"] = true;
    checkList["clusterSizeDist_exact"] = true;
    checkList["clusterSizeDist_time"] = true;
    checkList["interEventTimeDist_time"] =true;
    checkList["ageDist_time"] = true;
    checkList["deltaUpperBoundDist_time"] = true;
    checkList["orderParameterDist"] = true;


    //* Run data process
    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, logBinDelta, deletion);
    mBFW::data::process(checkList);
    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}