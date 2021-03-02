#include <iostream>
#include <chrono>

#include "parameters.hpp"
#include "data.hpp"

int main(int argc, char *argv[]){
    const double acceptanceThreshold=std::stod(argv[1]);
    const int networkSize=std::stoul(argv[2]);
    const bool deletion = false;
    constexpr double logBinDelta = 0.1;

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
    checkList["ageDist_op"] = false;
    checkList["ageDist_time"] = false;
    checkList["clusterSizeDist"] = false;
    checkList["deltaUpperBound_deltaAcceptance"] = false;
    checkList["deltaUpperBoundDist_op"] = false;
    checkList["deltaUpperBoundDist_time"] = false;
    checkList["deltaUpperBoundDist_tot"] = false;
    checkList["interEventTime"] = false;
    checkList["interEventTimeDist_op"] = false;
    checkList["interEventTimeDist_time"] = false;
    checkList["interEventTimeDist_tot"] = false;
    checkList["interEventTime_deltaUpperBound"] = false;
    checkList["meanClusterSize"] = false;
    checkList["orderParameter"] = false;
    checkList["orderParameterDist"] = false;
    checkList["orderParameterVariance"] = false;
    checkList["upperBound_deltaAcceptance"] = false;
    checkList["sampled_deltaUpperBound_interEventTime"] = true;
    checkList["sampled_upperBound_interEventTime"] = false;
    checkList["sampled_time_interEventTime"] = false;

    checkList["clusterSizeDist_exact"] = false;
    checkList["clusterSizeDist_time"] = false;
    checkList["meanClusterSize_trial"] = false;
    checkList["orderParameter_trial"] = false;
    checkList["orderParameterVariance_trial"] = false;

    //* Run data process
    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, logBinDelta, deletion);
    mBFW::data::process(checkList);
    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}