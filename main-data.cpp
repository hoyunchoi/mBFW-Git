#include <chrono>
#include <iostream>

#include "data.hpp"
#include "parameters.hpp"

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);
    const double acceptanceThreshold = std::stod(argv[1]);
    const int networkSize = std::stoul(argv[2]);
    mBFW::data::deletion = false;

    //* Check input network size and acceptance threshold
    if (mBFW::parameters::networkSizeList.find(networkSize) == mBFW::parameters::networkSizeList.end()) {
        std::cout << "WARNING: network size is not valid\n";
        return -1;
    }
    if (mBFW::parameters::acceptanceThresholdList.find(acceptanceThreshold) == mBFW::parameters::acceptanceThresholdList.end()) {
        std::cout << "WARNING: acceptance threshold is not valid\n";
        return -1;
    }

    //* Check list of each observables
    std::map<std::string, bool> checkList;
    checkList["ageDist_op"] = true;
    checkList["ageDist_time"] = true;
    checkList["clusterSizeDist"] = true;
    checkList["deltaUpperBound_deltaAcceptance"] = true;
    checkList["deltaUpperBoundDist_op"] = true;
    checkList["deltaUpperBoundDist_time"] = true;
    checkList["deltaUpperBoundDist_tot"] = false;
    checkList["interEventTime"] = true;
    checkList["interEventTimeDist_op"] = true;
    checkList["interEventTimeDist_time"] = true;
    checkList["interEventTimeDist_tot"] = false;
    checkList["interEventTime_deltaUpperBound"] = false;
    checkList["meanClusterSize"] = true;
    checkList["orderParameter"] = true;
    checkList["orderParameterDist"] = true;
    checkList["orderParameterVariance"] = true;
    checkList["upperBound_deltaAcceptance"] = true;
    checkList["sampled_deltaUpperBound_interEventTime"] = false;
    checkList["sampled_upperBound_interEventTime"] = false;
    checkList["sampled_time_interEventTime"] = false;
    checkList["clusterSizeDist_exact"] = true;
    checkList["clusterSizeDist_time"] = true;
    checkList["meanClusterSize_trial"] = true;
    checkList["orderParameter_trial"] = true;
    checkList["orderParameterVariance_trial"] = true;

    //* Run
    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold);
    mBFW::data::run(checkList);
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    printf("%.6f second to process data\n", sec.count());

    return 0;
}