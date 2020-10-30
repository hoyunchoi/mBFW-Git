#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <filesystem>

#include "../library-Git/linearAlgebra.hpp"
#include "../library-Git/Networks.hpp"
#include "../library-Git/CSV.hpp"
#include "../library-Git/pcg_random.hpp"

#include "parameters.hpp"
#include "mBFW.hpp"

namespace mBFW::generate{
    //*--------------------------------------------Declaration of variables---------------------------------------------------------
    //* Declaration of variables only used at mBFW::generate namespace
    int randomEngineSeed;
    double precision;
    double degenerated;

    //* Declaration of Random Engine
    pcg32 randomEngine;
    std::uniform_int_distribution<int> nodeDistribution;

    //* Declaration of observables
    //! orderParameter[time]: Average value of order parameter at specific time
    std::vector<double> orderParameter;
    std::vector<double> orderParameter_trial;

    //! secondMoment[time]: Average value of second moment at specific time
    std::vector<double> secondMoment;
    std::vector<double> secondMoment_trial;

    //! meanClusterSize[time]: Average value of normalized mean cluster size at specific time
    std::vector<double> meanClusterSize;
    std::vector<double> meanClusterSize_trial;

    //! clusterSizeDist[op]: Average distribution of cluster size when order parameter passes op
    //! clusterSizeDist_exact[op]: Average distribution of cluster size when order parameter is exactly op
    //! clusterSizeDist_time[time]: Average distribution of cluster size when time at specific time
    std::map<double, std::vector<long long>> clusterSizeDist;
    std::map<double, std::vector<long long>> clusterSizeDist_exact;
    std::map<double, std::vector<long long>> clusterSizeDist_time;


    //*-------------------------------------------Set Parameters for one run------------------------------------------------------
    void setParameters(const int& t_networkSize, const int& t_ensembleSize, const double& t_acceptanceThreshold, const int& t_coreNum, const int& t_randomEngineSeed){
        //* Input variables
        networkSize = t_networkSize;
        ensembleSize = t_ensembleSize;
        coreNum = t_coreNum;
        acceptanceThreshold = t_acceptanceThreshold;
        randomEngineSeed = t_randomEngineSeed;
        orderParameter_clusterSizeDist = mBFW::parameters::set_orderParameter_clusterSizeDist(t_networkSize, t_ensembleSize);
        time_clusterSizeDist = mBFW::parameters::set_time_clusterSizeDist(t_networkSize, t_ensembleSize);
        precision = 1e4;
        if (t_networkSize < precision){
            precision = t_networkSize;
        }

        //* Initialize Random Engine
        randomEngineSeed == -1 ? randomEngine.seed((std::random_device())()) : randomEngine.seed(randomEngineSeed);
        nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, networkSize-1));

        //* Resize observables
        const int maxTime = t_networkSize;
        const int maxTrialTime = std::floor(maxTime/t_acceptanceThreshold);
        orderParameter.assign(maxTime, 0.0);
        orderParameter_trial.assign(maxTrialTime, 0.0);
        secondMoment.assign(maxTime, 0.0);
        secondMoment_trial.assign(maxTrialTime, 0.0);
        meanClusterSize.assign(maxTime, 0.0);
        meanClusterSize_trial.assign(maxTrialTime, 0.0);
        for (const double& op : orderParameter_clusterSizeDist){
            clusterSizeDist[op].assign(t_networkSize, 0);
            clusterSizeDist_exact[op].assign(t_networkSize, 0);
        }
        for (const double& t : time_clusterSizeDist){
            clusterSizeDist_time[t].assign(t_networkSize, 0);
        }
    } //* End of function mBFW::generate::setParameters

    //*-------------------------------------------Run mBFW model ------------------------------------------------------
    void run(){
        for (int ensemble=0; ensemble<ensembleSize; ++ensemble){
            //* Default values for one ensemble
            NZ_Network model(networkSize);
            Node node1, node2, root1, root2;
            int size1, size2;
            int time = 0;
            int trialTime = 0;
            int upperBound = 2;
            int eventTime = 0;
            bool choosingNewNode = true;
            std::set<double> findingClusterSizeDist;
            for (const double& op : orderParameter_clusterSizeDist){
                findingClusterSizeDist.insert(op);
            }
            std::set<double> newFindingClusterSizeDist = findingClusterSizeDist;

            //* initial condition
            orderParameter[0] += 1.0/networkSize;
            orderParameter_trial[0] += 1.0/networkSize;
            secondMoment[0] += std::pow(1.0/networkSize, 2.0);
            secondMoment_trial[0] += std::pow(1.0/networkSize, 2.0);
            meanClusterSize[0] += 1.0;
            meanClusterSize_trial[0] += 1.0;

            //* Do rBFW algorithm until all clusters merge to one
            while (model.getMaximumClusterSize() < networkSize){
                //* Find new nodes
                if (choosingNewNode){
                    //* Randomly choose new nodes
                    do {
                        node1 = nodeDistribution(randomEngine);
                        node2 = nodeDistribution(randomEngine);
                        root1 = model.getRoot(node1);
                        root2 = model.getRoot(node2);
                    } while(root1 == root2);

                    //* choose two clusters of each node
                    size1 = model.getClusterSize(root1);
                    size2 = model.getClusterSize(root2);
                }

                //* Merge two clusters, update time
                if (size1+size2 <= upperBound){
                    model.merge(root1, root2);
                    ++time;
                    ++trialTime;
                    choosingNewNode = true;
                    const int currentMaximumClusterSize = model.getMaximumClusterSize();
                    const double currentOrderParameter = (double)currentMaximumClusterSize/networkSize;

                    //! Order Parameter
                    orderParameter[time] += currentOrderParameter;
                    orderParameter_trial[trialTime] += currentOrderParameter;

                    //! Second Moment
                    secondMoment[time] += std::pow(currentOrderParameter, 2.0);
                    secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);

                    //! Mean Cluster Size
                    meanClusterSize[time] += model.getMeanClusterSize();
                    meanClusterSize_trial[trialTime] += model.getMeanClusterSize();

                    //! Cluster Size Distribution time
                    auto it = std::find(time_clusterSizeDist.begin(), time_clusterSizeDist.end(), (double)time/networkSize);
                    if (it != time_clusterSizeDist.end()){
                        const std::map<int, int> sortedCluster = model.getSortedCluster();
                        for (auto it2 = sortedCluster.begin(); it2 != sortedCluster.end(); ++it2){
                            clusterSizeDist_time[*it][it2->first] += it2->second;
                        }
                    }

                    // //* Order Parameter of network is changed
                    if (model.getDeltaMaximumClusterSize()){
                        //! Cluster Size Distribution
                        findingClusterSizeDist = newFindingClusterSizeDist;
                        for (const double& op : findingClusterSizeDist){
                            if (op < currentOrderParameter){
                                const std::map<int, int> sortedCluster = model.getSortedCluster();
                                for (auto it2 = sortedCluster.begin(); it2 != sortedCluster.end(); ++it2){
                                    clusterSizeDist[op][it2->first] += it2->second;
                                }
                                newFindingClusterSizeDist.erase(op);
                            }
                        }

                        //! Cluster Size Distribution Exact
                        const double roundedOrderParameter = round(currentOrderParameter*precision)/precision;
                        auto it = std::find(orderParameter_clusterSizeDist.begin(), orderParameter_clusterSizeDist.end(), roundedOrderParameter);
                        if (it != orderParameter_clusterSizeDist.end()){
                            const std::map<int, int> sortedCluster = model.getSortedCluster();
                            for (auto it2 = sortedCluster.begin(); it2 != sortedCluster.end(); ++it2){
                                clusterSizeDist_exact[*it][it2->first] += it2 -> second;
                            }
                        }
                    }//* End of order parameter update
                }

                //* Upper Bound change
                else if ((double)time/trialTime <= acceptanceThreshold){
                    upperBound = size1+size2;
                    choosingNewNode=false;
                }//* End of Upper Bound change

                //* Choosed link rejected
                else{
                    ++trialTime;
                    choosingNewNode=true;
                    const double currentOrderParameter = (double)model.getMaximumClusterSize()/networkSize;

                    //! Trial Order Parameter
                    orderParameter_trial[trialTime] += currentOrderParameter;

                    //! Trial Second Moment
                    secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);

                    //! Trial Mean Cluster Size
                    meanClusterSize_trial[trialTime] += model.getMeanClusterSize();
                }//* End of rejecting link
                //* End of one step
            }//* End of network growing (one ensemble)
        } //* End of every ensembles
    } //* End of function mBFW::generate::run

    //*-------------------------------------------Save calculated variables------------------------------------------------------
    void save(){
        using namespace linearAlgebra;
        namespace fs = std::filesystem;

        //! Average and Save Order Parameter
        orderParameter /= ensembleSize;
        const std::string orderParameterPath = rootPath + "orderParameter/";
        if (!fs::exists(orderParameterPath)){
            fs::create_directories(orderParameterPath);
        }
        CSV::write(orderParameterPath + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameter);

        //! Average and Save Trial Order Parameter
        orderParameter_trial /= ensembleSize;
        const std::string orderParameter_trialPath = rootPath + "orderParameter_trial/";
        if (!fs::exists(orderParameter_trialPath)){
            fs::create_directories(orderParameter_trialPath);
        }
        CSV::write(orderParameter_trialPath + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameter_trial);

        //! Average and Save Order Parameter Variance
        secondMoment /= ensembleSize;
        const std::vector<double> orderParameterVariance = secondMoment - elementPow(orderParameter, 2.0);
        const std::string orderParameterVariancePath = rootPath + "orderParameterVariance/";
        if (!fs::exists(orderParameterVariancePath)){
            fs::create_directories(orderParameterVariancePath);
        }
        CSV::write(orderParameterVariancePath + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameterVariance);

        //! Average and Save Order Parameter Trial Variance
        secondMoment_trial /= ensembleSize;
        const std::vector<double> orderParameterVariance_trial = secondMoment_trial - elementPow(orderParameter_trial, 2.0);
        const std::string orderParameterVariance_trialPath = rootPath + "orderParameterVariance_trial/";
        if (!fs::exists(orderParameterVariance_trialPath)){
            fs::create_directories(orderParameterVariance_trialPath);
        }
        CSV::write(orderParameterVariance_trialPath + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameterVariance_trial);

        //! Average and Save Mean Cluster Size
        meanClusterSize /= ensembleSize;
        const std::string meanClusterSizePath = rootPath + "meanClusterSize/";
        if (!fs::exists(meanClusterSizePath)){
            fs::create_directories(meanClusterSizePath);
        }
        CSV::write(meanClusterSizePath + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum), meanClusterSize);

        //! Average and Save Trial Mean Cluster Size
        meanClusterSize_trial /= ensembleSize;
        const std::string meanClusterSize_trialPath = rootPath + "meanClusterSize_trial/";
        if (!fs::exists(meanClusterSize_trialPath)){
            fs::create_directories(meanClusterSize_trialPath);
        }
        CSV::write(meanClusterSize_trialPath + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum), meanClusterSize_trial);

        //! Cluster Size Distribution
        const std::string clusterSizeDistPath = rootPath + "clusterSizeDist/";
        if (!fs::exists(clusterSizeDistPath)){
            fs::create_directories(clusterSizeDistPath);
        }
        for (const double& op : orderParameter_clusterSizeDist){
            std::map<int, double> trimmed;
            const long long total = std::accumulate(clusterSizeDist[op].begin(),clusterSizeDist[op].end(), 0);
            for (int cs=0; cs<networkSize; ++cs){
                if (clusterSizeDist[op][cs]){
                    trimmed[cs] = (double)clusterSizeDist[op][cs]/total;
                }
            }
            CSV::write(clusterSizeDistPath + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSize, op, coreNum), trimmed);
        }

        //! Cluster Size Distribution Exact
        const std::string clusterSizeDist_exactPath = rootPath + "clusterSizeDist_exact/";
        if (!fs::exists(clusterSizeDist_exactPath)){
            fs::create_directories(clusterSizeDist_exactPath);
        }
        for (const double& op : orderParameter_clusterSizeDist){
            std::map<int, double> trimmed;
            const long long total = std::accumulate(clusterSizeDist_exact[op].begin(),clusterSizeDist_exact[op].end(), 0);
            for (int cs=0; cs<networkSize; ++cs){
                if (clusterSizeDist_exact[op][cs]){
                    trimmed[cs] = (double)clusterSizeDist_exact[op][cs]/total;
                }
            }
            CSV::write(clusterSizeDist_exactPath + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSize, op, coreNum), trimmed);
        }

        //! Cluster Size Distribution Time
        const std::string clusterSizeDist_timePath = rootPath + "clusterSizeDist_time/";
        if (!fs::exists(clusterSizeDist_timePath)){
            fs::create_directories(clusterSizeDist_timePath);
        }
        for (const double& t : time_clusterSizeDist){
            std::map<int, double> trimmed;
            const long long total = std::accumulate(clusterSizeDist_time[t].begin(), clusterSizeDist_time[t].end(), 0);
            for (int cs=0; cs<networkSize; ++cs){
                if (clusterSizeDist_time[t][cs]){
                    trimmed[cs] = (double)clusterSizeDist_time[t][cs]/total;
                }
            }
            CSV::write(clusterSizeDist_timePath + filename_time(networkSize, acceptanceThreshold, ensembleSize, t, coreNum), trimmed);
        }
    }//* End of function mBFW::generate::save
}//* End of namespace mBFW::generate