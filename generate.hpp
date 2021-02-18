#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <random>
#include <algorithm>
#include <functional>
#include <filesystem>

#include "../library/linearAlgebra.hpp"
#include "../library/CSV.hpp"
#include "../library/pcg_random.hpp"

#include "NZ_Network.hpp"
#include "parameters.hpp"
#include "fileName.hpp"

namespace mBFW::generate{
    const std::string baseDirectory = "../data/mBFW/";

    //*--------------------------------------------Declaration of variables---------------------------------------------------------
    //* Declaration of variables used at mBFW::generate namespace
    int networkSize;
    double acceptanceThreshold;
    int ensembleSize;
    int coreNum;
    int randomEngineSeed;
    double precision;
    int maxTime;
    int maxTrialTime;
    const std::vector<std::string> states = {"before", "during", "after"};
    double t_a, m_a, t_c, m_c;

    //* Declaration of Random Engine
    pcg32 randomEngine;
    std::uniform_int_distribution<int> nodeDistribution;

    //* Declaration of observables
    //! orderParameter[time]: Average value of order parameter at specific time
    std::vector<double> orderParameter;
    // std::vector<double> orderParameter_trial;

    //! secondMoment[time]: Average value of second moment at specific time
    std::vector<double> secondMoment;
    // std::vector<double> secondMoment_trial;

    //! meanClusterSize[time]: Average value of normalized mean cluster size at specific time
    std::vector<double> meanClusterSize;
    // std::vector<double> meanClusterSize_trial;

    //! interEventTime[time]: Average value of inter event time at specific time
    std::vector<double> interEventTime;
    std::vector<int> sampled_interEventTime;

    //! clusterSizeDist[op]: Average distribution of cluster size when order parameter passes op
    //! clusterSizeDist_exact[op]: Average distribution of cluster size when order parameter is exactly op
    //! clusterSizeDist_time[time]: Average distribution of cluster size when time at specific time
    std::set<double> orderParameter_clusterSizeDist;
    // std::set<double> time_clusterSizeDist;
    std::map<double, std::vector<long long>> clusterSizeDist;
    // std::map<double, std::vector<long long>> clusterSizeDist_exact;
    // std::map<double, std::vector<long long>> clusterSizeDist_time;

    //! orderParameterDist[time] : Distribution of order parameter at specific time
    // std::set<double> time_orderParameterDist;
    // std::map<double, std::vector<int>> orderParameterDist;

    //! interEventTimeDist_op(time)["before"] : Average distribution of inter event time before discontinuous jump by order paramter (time)
    //! interEventTimeDist_op(time)["during"] : Average distribution of inter event time during discontinuous jump by order paramter (time)
    //! interEventTimeDist_op(time)["after"] : Average distribution of inter event time after discontinuous jump by order paramter (time)
    std::map<std::string, std::vector<int>> interEventTimeDist_op;
    std::map<std::string, std::vector<int>> interEventTimeDist_time;

    //! ageDist_op(time)["state"] : Average distribution of age before(during, after) discontinuous jump by order paramter (time)
    std::map<std::string, std::vector<long long>> ageDist_op;
    std::map<std::string, std::vector<long long>> ageDist_time;

    //! deltaUpperBoundDist_op(time)["state"] : Average distribution of deltaK before(during, after) discontinuous jump by order paramter (time)
    std::map<std::string, std::vector<int>> deltaUpperBoundDist_op;
    std::map<std::string, std::vector<int>> deltaUpperBoundDist_time;

    //! upperBound_deltaAcceptance[state][K] : Average delta acceptance at specific upper bound K before(during, after) discontinuous jump by order parameter
    std::map<std::string, std::vector<double>> upperBound_deltaAcceptance;
    std::map<std::string, std::vector<int>> sampled_upperBound_deltaAcceptance;

    //! deltaUpperBound_deltaAcceptance[state][deltaK] : Average delta acceptance at specific delta upper bound deltaK before(during, after) discontinuous jump by order parameter
    std::map<std::string, std::vector<double>> deltaUpperBound_deltaAcceptance;
    std::map<std::string, std::vector<int>> sampled_deltaUpperBound_deltaAcceptance;

    //*-------------------------------------------Set Parameters for one run------------------------------------------------------
    void setParameters(const int& t_networkSize, const int& t_ensembleSize, const double& t_acceptanceThreshold, const int& t_coreNum, const int& t_randomEngineSeed){
        //* Input variables
        networkSize = t_networkSize;
        ensembleSize = t_ensembleSize;
        coreNum = t_coreNum;
        acceptanceThreshold = t_acceptanceThreshold;
        randomEngineSeed = t_randomEngineSeed;
        orderParameter_clusterSizeDist = mBFW::parameters::set_orderParameter_clusterSizeDist(t_networkSize, t_acceptanceThreshold);
        // time_clusterSizeDist = mBFW::parameters::set_time_clusterSizeDist(t_networkSize, t_acceptanceThreshold);
        // time_orderParameterDist = mBFW::parameters::set_time_orderParameterDist(t_networkSize, t_acceptanceThreshold);
        std::tie(t_a, m_a, t_c, m_c) = mBFW::parameters::set_points(t_networkSize, t_acceptanceThreshold);
        t_a *= t_networkSize; t_c *= t_networkSize;
        precision = 1e4;
        if (t_networkSize < precision){
            precision = t_networkSize;
        }

        //* Initialize Random Engine
        randomEngineSeed == -1 ? randomEngine.seed((std::random_device())()) : randomEngine.seed(randomEngineSeed);
        nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, networkSize-1));

        //* Resize observables
        maxTime = t_networkSize;
        maxTrialTime = std::floor(maxTime/t_acceptanceThreshold);
        orderParameter.assign(maxTime, 0.0);
        secondMoment.assign(maxTime, 0.0);
        meanClusterSize.assign(maxTime, 0.0);
        // orderParameter_trial.assign(maxTrialTime, 0.0);
        // secondMoment_trial.assign(maxTrialTime, 0.0);
        // meanClusterSize_trial.assign(maxTrialTime, 0.0);
        interEventTime.assign(maxTime, 0.0);
        sampled_interEventTime.assign(maxTime, 0);
        for (const double& op : orderParameter_clusterSizeDist){
            clusterSizeDist[op].assign(t_networkSize, 0);
            // clusterSizeDist_exact[op].assign(t_networkSize, 0);
        }
        // for (const double& t : time_clusterSizeDist){
        //     clusterSizeDist_time[t].assign(t_networkSize, 0);
        // }
        // for (const double& t : time_orderParameterDist){
        //     orderParameterDist[t].assign(t_networkSize, 0);
        // }
        for (const std::string& state : states){
            interEventTimeDist_op[state].assign(t_networkSize, 0);
            interEventTimeDist_time[state].assign(t_networkSize, 0);
            ageDist_op[state].assign(t_networkSize, 0);
            ageDist_time[state].assign(t_networkSize, 0);
            deltaUpperBoundDist_op[state].assign(t_networkSize, 0);
            deltaUpperBoundDist_time[state].assign(t_networkSize, 0);
            upperBound_deltaAcceptance[state].assign(t_networkSize, 0.0);
            sampled_upperBound_deltaAcceptance[state].assign(t_networkSize, 0);
            deltaUpperBound_deltaAcceptance[state].assign(t_networkSize, 0.0);
            sampled_deltaUpperBound_deltaAcceptance[state].assign(t_networkSize, 0);
        }
    } //* End of function mBFW::generate::setParameters

    //*-------------------------------------------Run mBFW model ------------------------------------------------------
    void run(){
        for (int ensemble=0; ensemble<ensembleSize; ++ensemble){
            //* Default values for one ensemble
            NZ_Network model(networkSize);
            int root1, root2;
            std::string currentState_time = "before";
            std::string currentState_op = "before";
            int size1, size2;
            int time = 0;
            int trialTime = 0;
            int upperBound = 2;
            int eventTime = 0;
            double maxDeltaAcceptance = 0.0;
            bool choosingNewNode = true;
            std::set<double> findingClusterSizeDist;
            for (const double& op : orderParameter_clusterSizeDist){
                findingClusterSizeDist.insert(op);
            }
            std::set<double> newFindingClusterSizeDist = findingClusterSizeDist;

            //* Set initial conditions
            orderParameter[0] += 1.0/networkSize;
            secondMoment[0] += std::pow(1.0/networkSize, 2.0);
            meanClusterSize[0] += 1.0;
            // orderParameter_trial[0] += 1.0/networkSize;
            // secondMoment_trial[0] += std::pow(1.0/networkSize, 2.0);
            // meanClusterSize_trial[0] += 1.0;

            //* Do mBFW algorithm until all clusters merge to one
            while (model.getMaximumClusterSize() < networkSize){
                //* Find new nodes
                if (choosingNewNode){
                    //* Randomly choose new nodes
                    do {
                        root1 = model.getRoot(nodeDistribution(randomEngine));
                        root2 = model.getRoot(nodeDistribution(randomEngine));
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

                    //* Check the state distinguished by time
                    if (time < t_a){
                        currentState_time = "before";
                    }
                    else if (time < t_c){
                        currentState_time = "during";
                    }
                    else{
                        currentState_time = "after";
                    }

                    //* Update max delta acceptance
                    maxDeltaAcceptance = std::max((double)time/trialTime - acceptanceThreshold, maxDeltaAcceptance);

                    //! Order Parameter
                    {
                        orderParameter[time] += currentOrderParameter;
                        // orderParameter_trial[trialTime] += currentOrderParameter;
                    }

                    //! Second Moment
                    {
                        secondMoment[time] += std::pow(currentOrderParameter, 2.0);
                        // secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);
                    }

                    //! Mean Cluster Size
                    {
                        meanClusterSize[time] += model.getMeanClusterSize();
                        // meanClusterSize_trial[trialTime] += model.getMeanClusterSize();
                    }

                    //! Cluster Size Distribution time
                    {
                        // auto it = std::find(time_clusterSizeDist.begin(), time_clusterSizeDist.end(), (double)time/networkSize);
                        // if (it != time_clusterSizeDist.end()){
                        //     const std::map<int, int> sortedCluster = model.getSortedCluster();
                        //     for (auto it2 = sortedCluster.begin(); it2 != sortedCluster.end(); ++it2){
                        //         clusterSizeDist_time[*it][it2->first] += it2->second;
                        //     }
                        // }
                    }

                    //! Order Parameter Distribution
                    {
                        // auto it = std::find(time_orderParameterDist.begin(), time_orderParameterDist.end(), (double)time/networkSize);
                        // if (it != time_orderParameterDist.end()){
                        //     ++orderParameterDist[*it][currentMaximumClusterSize];
                        // }
                    }

                    //! Age Distribution time(op)
                    {
                        const std::vector<std::pair<int, int>> changedAge = model.getChangedAge();
                        for (const auto& age : changedAge){
                            ageDist_time[currentState_time][age.first] += age.second;
                            ageDist_op[currentState_op][age.first] += age.second;
                        }
                    }

                    //* Order Parameter of network is changed <=> upper bound is changed
                    if (model.getDeltaMaximumClusterSize() && currentMaximumClusterSize > 2){
                        const int deltaMaximumClusterSize = model.getDeltaMaximumClusterSize();
                        //* Check the state distinguished by order parameter
                        if (currentOrderParameter < m_a){
                            currentState_op = "before";
                        }
                        else if (currentOrderParameter < m_c){
                            currentState_op = "during";
                        }
                        else{
                            currentState_op = "after";
                        }

                        //! Cluster Size Distribution
                        {
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
                        }

                        //! Cluster Size Distribution Exact
                        {
                            // const double roundedOrderParameter = round(currentOrderParameter*precision)/precision;
                            // auto it = std::find(orderParameter_clusterSizeDist.begin(), orderParameter_clusterSizeDist.end(), roundedOrderParameter);
                            // if (it != orderParameter_clusterSizeDist.end()){
                            //     const std::map<int, int> sortedCluster = model.getSortedCluster();
                            //     for (auto it2 = sortedCluster.begin(); it2 != sortedCluster.end(); ++it2){
                            //         clusterSizeDist_exact[*it][it2->first] += it2 -> second;
                            //     }
                            // }
                        }

                        //! Inter Event Time
                        {
                            interEventTime[time] += (double)time-eventTime;
                            ++sampled_interEventTime[time];
                        }

                        //! Inter Event Time Distribution time(op)
                        {
                            ++interEventTimeDist_time[currentState_time][time-eventTime];
                            ++interEventTimeDist_op[currentState_op][time-eventTime];

                        }

                        //! Delta Upper Bound Distribution time(op)
                        {
                            ++deltaUpperBoundDist_time[currentState_time][deltaMaximumClusterSize];
                            ++deltaUpperBoundDist_op[currentState_op][deltaMaximumClusterSize];
                        }

                        //! (Delta)Upper Bound_Delta Acceptance
                        {
                            upperBound_deltaAcceptance[currentState_op][upperBound] += maxDeltaAcceptance;
                            ++sampled_upperBound_deltaAcceptance[currentState_op][upperBound];
                            deltaUpperBound_deltaAcceptance[currentState_op][deltaMaximumClusterSize] += maxDeltaAcceptance;
                            ++sampled_deltaUpperBound_deltaAcceptance[currentState_op][upperBound];
                        }

                        //* Initialize variable for new period
                        eventTime = time;
                        maxDeltaAcceptance = 0;
                    }//* End of order parameter update, End of k-period
                }

                //* Upper Bound change
                else if ((double)time/trialTime <= acceptanceThreshold){
                    upperBound = size1+size2;
                    choosingNewNode=false;
                }//* End of Upper Bound change

                //* Chosen link rejected
                else{
                    ++trialTime;
                    choosingNewNode=true;
                    // const double currentOrderParameter = (double)model.getMaximumClusterSize()/networkSize;

                    //! Trial Order Parameter
                    // orderParameter_trial[trialTime] += currentOrderParameter;

                    //! Trial Second Moment
                    // secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);

                    //! Trial Mean Cluster Size
                    // meanClusterSize_trial[trialTime] += model.getMeanClusterSize();
                }//* End of rejecting link
                //* End of one step
            }//* End of network growing (one ensemble)
        } //* End of every ensembles
    } //* End of function mBFW::generate::run

    //*-------------------------------------------Save calculated variables------------------------------------------------------
    void save(){
        using namespace linearAlgebra;
        namespace fs = std::filesystem;

        //! Order Parameter
        {
            orderParameter /= ensembleSize;
            const std::string directory = baseDirectory + "orderParameter/";
            if (!fs::exists(directory)){
                fs::create_directories(directory);
            }
            CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameter);
        }

        //! Trial Order Parameter
        {
            // orderParameter_trial /= ensembleSize;
            // const std::string directory = baseDirectory + "orderParameter_trial/";
            // if (!fs::exists(directory)){
            //     fs::create_directories(directory);
            // }
            // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameter_trial);
        }

        //! Order Parameter Variance
        {
            secondMoment /= ensembleSize;
            std::vector<double> orderParameterVariance(maxTime, 0.0);
            for (int i=0; i<maxTime; ++i){
                const double var = secondMoment[i] - pow(orderParameter[i], 2.0);
                var < 0 ? orderParameterVariance[i] = 0.0 :orderParameterVariance[i] = std::sqrt(var) * networkSize;
            }
            const std::string directory = baseDirectory + "orderParameterVariance/";
            if (!fs::exists(directory)){
                fs::create_directories(directory);
            }
            CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameterVariance);
        }

        //! Trial Order Parameter Variance
        {
            // secondMoment_trial /= ensembleSize;
            // std::vector<double> orderParameterVariance_trial(maxTrialTime, 0.0);
            // for (int i=0; i<maxTrialTime; ++i){
            //     const double var = secondMoment_trial[i] - pow(orderParameter_trial[i], 2.0);
            //     var<0 ? orderParameterVariance_trial[i] = 0.0 : orderParameterVariance_trial[i] = std::sqrt(var) * networkSize;
            // }
            // const std::string directory = baseDirectory + "orderParameterVariance_trial/";
            // if (!fs::exists(directory)){
            //     fs::create_directories(directory);
            // }
            // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameterVariance_trial);
        }

        //! Mean Cluster Size
        {
            meanClusterSize /= ensembleSize;
            const std::string directory = baseDirectory + "meanClusterSize/";
            if (!fs::exists(directory)){
                fs::create_directories(directory);
            }
            CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), meanClusterSize);
        }

        //! Trial Mean Cluster Size
        {
            // meanClusterSize_trial /= ensembleSize;
            // const std::string directory = baseDirectory + "meanClusterSize_trial/";
            // if (!fs::exists(directory)){
            //     fs::create_directories(directory);
            // }
            // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), meanClusterSize_trial);
        }

        //! Order Parameter Distribution
        {
            // const std::string directory = baseDirectory + "orderParameterDist/";
            // if(!fs::exists(directory)){
            //     fs::create_directories(directory);
            // }
            // for (const double& t : time_orderParameterDist){
            //     std::map<double, double> trimmed;
            //     for (int op=0; op<networkSize; ++op){
            //         if (orderParameterDist[t][op]){
            //             trimmed[(double)op/networkSize] = (double)orderParameterDist[t][op] / ensembleSize;
            //         }
            //     }
            //     CSV::write(directory + fileName::NGET(networkSize, acceptanceThreshold, ensembleSize, t, coreNum), trimmed);
            // }
        }

        //! Cluster Size Distribution
        {
            const std::string directory = baseDirectory + "clusterSizeDist/";
            if (!fs::exists(directory)){
                fs::create_directories(directory);
            }
            for (const double& op : orderParameter_clusterSizeDist){
                std::map<int, double> trimmed;
                const long long tot = std::accumulate(clusterSizeDist[op].begin(),clusterSizeDist[op].end(), 0);
                for (int cs=0; cs<networkSize; ++cs){
                    if (clusterSizeDist[op][cs]){
                        trimmed[cs] = (double)clusterSizeDist[op][cs] / tot;
                    }
                }
                CSV::write(directory + fileName::NGEOP(networkSize, acceptanceThreshold, ensembleSize, op, coreNum), trimmed);
            }
        }

        //! Cluster Size Distribution Exact
        {
            // const std::string directory = baseDirectory + "clusterSizeDist_exact/";
            // if (!fs::exists(directory)){
            //     fs::create_directories(directory);
            // }
            // for (const double& op : orderParameter_clusterSizeDist){
            //     std::map<int, double> trimmed;
            //     const long long tot = std::accumulate(clusterSizeDist_exact[op].begin(),clusterSizeDist_exact[op].end(), 0);
            //     for (int cs=0; cs<networkSize; ++cs){
            //         if (clusterSizeDist_exact[op][cs]){
            //             trimmed[cs] = (double)clusterSizeDist_exact[op][cs]/tot;
            //         }
            //     }
            //     CSV::write(directory + fileName::NGEOP(networkSize, acceptanceThreshold, ensembleSize, op, coreNum), trimmed);
            // }
        }

        //! Cluster Size Distribution Time
        {
            // const std::string directory = baseDirectory + "clusterSizeDist_time/";
            // if (!fs::exists(directory)){
            //     fs::create_directories(directory);
            // }
            // for (const double& t : time_clusterSizeDist){
            //     std::map<int, double> trimmed;
            //     const long long tot = std::accumulate(clusterSizeDist_time[t].begin(), clusterSizeDist_time[t].end(), 0);
            //     for (int cs=0; cs<networkSize; ++cs){
            //         if (clusterSizeDist_time[t][cs]){
            //             trimmed[cs] = (double)clusterSizeDist_time[t][cs]/tot;
            //         }
            //     }
            //     CSV::write(directory + fileName::NGET(networkSize, acceptanceThreshold, ensembleSize, t, coreNum), trimmed);
            // }
        }

        //! Inter Event Time
        {
            const std::string directory = baseDirectory + "interEventTime/";
            if (!fs::exists(directory)){
                fs::create_directories(directory);
            }
            for (int t=0; t<maxTime; ++t){
                if (sampled_interEventTime[t]){
                    interEventTime[t] /= (double)sampled_interEventTime[t];
                }
            }
            CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), interEventTime);
        }

        //! Inter Event Time Distribution time(op)
        {
            for (const std::string& state : states){
                const std::string directory_time = baseDirectory + "interEventTimeDist_time/" + state + "/";
                if (!fs::exists(directory_time)){
                    fs::create_directories(directory_time);
                }
                std::map<int, double> trimmed_time;
                const double tot_time = std::accumulate(interEventTimeDist_time[state].begin(), interEventTimeDist_time[state].end(), 0.0);
                for (int iet=0; iet<networkSize; ++iet){
                    if (interEventTimeDist_time[state][iet]){
                        trimmed_time[iet] = interEventTimeDist_time[state][iet] / tot_time;
                    }
                }
                CSV::write(directory_time + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_time);

                const std::string directory_op = baseDirectory + "interEventTimeDist_op/" + state + "/";
                if (!fs::exists(directory_op)){
                    fs::create_directories(directory_op);
                }
                std::map<int, double> trimmed_op;
                const double tot_op = std::accumulate(interEventTimeDist_op[state].begin(), interEventTimeDist_op[state].end(), 0.0);
                for (int iet=0; iet<networkSize; ++iet){
                    if (interEventTimeDist_op[state][iet]){
                        trimmed_op[iet] = interEventTimeDist_op[state][iet] / tot_op;
                    }
                }
                CSV::write(directory_op + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_op);
            }
        }

        //! Age Distribution time(op)
        {
            for (const auto& state : states){
                const std::string directory_time = baseDirectory + "ageDist_time/" + state + "/";
                if (!fs::exists(directory_time)){
                    fs::create_directories(directory_time);
                }
                std::vector<double> temp_time(networkSize, 0.0);
                const double tot_time = std::accumulate(ageDist_time[state].begin(), ageDist_time[state].end(), 0.0);
                for (int age=0; age<networkSize; ++age){
                    temp_time[age] = ageDist_time[state][age] / tot_time;
                }
                print(temp_time);
                CSV::write(directory_time + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), temp_time);

                const std::string directory_op = baseDirectory + "ageDist_op/" + state + "/";
                if (!fs::exists(directory_op)){
                    fs::create_directories(directory_op);
                }
                std::vector<double> temp_op(networkSize, 0.0);
                const double tot_op = std::accumulate(ageDist_op[state].begin(), ageDist_op[state].end(), 0.0);
                for (int age=0; age<networkSize; ++age){
                    temp_op[age] = ageDist_time[state][age] / tot_time;
                }
                CSV::write(directory_op + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), temp_op);
            }
        }

        //! Delta Upper Bound Distribution_time(op)
        {
            for (const auto& state : states){
                const std::string directory_time = baseDirectory + "deltaUpperBoundDist_time/" + state + "/";
                if (!fs::exists(directory_time)){
                    fs::create_directories(directory_time);
                }
                std::map<int, double> trimmed_time;
                const double tot_time = std::accumulate(deltaUpperBoundDist_time[state].begin(), deltaUpperBoundDist_time[state].end(), 0.0);
                for (int deltaK = 0; deltaK < networkSize; ++deltaK){
                    if (deltaUpperBoundDist_time[state][deltaK]){
                        trimmed_time[deltaK] = (double)deltaUpperBoundDist_time[state][deltaK] / tot_time;
                    }
                }
                CSV::write(directory_time + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_time);

                const std::string directory_op = baseDirectory + "deltaUpperBoundDist_op/" + state + "/";
                if (!fs::exists(directory_op)){
                    fs::create_directories(directory_op);
                }
                std::map<int, double> trimmed_op;
                const double tot_op = std::accumulate(deltaUpperBoundDist_op[state].begin(), deltaUpperBoundDist_op[state].end(), 0.0);
                for (int deltaK=0; deltaK<networkSize; ++deltaK){
                    if (deltaUpperBoundDist_op[state][deltaK]){
                        trimmed_op[deltaK] = (double)deltaUpperBoundDist_op[state][deltaK] / tot_op;
                    }
                }
                CSV::write(directory_op + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_op);
            }
        }

        //! (Delta) Upper Bound vs Delta Acceptance
        {
            for (const auto& state : states){
                const std::string directory_k = baseDirectory + "upperBound_deltaAcceptance/" + state + "/";
                if (!fs::exists(directory_k)){
                    fs::create_directories(directory_k);
                }
                std::map<int, double> trimmed_k;
                for (int k=0; k<networkSize; ++k){
                    if (sampled_upperBound_deltaAcceptance[state][k] && upperBound_deltaAcceptance[state][k]){
                        trimmed_k[k] = upperBound_deltaAcceptance[state][k] / sampled_upperBound_deltaAcceptance[state][k];
                    }
                }
                CSV::write(directory_k + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_k);


                const std::string directory_deltaK = baseDirectory + "deltaUpperBound_deltaAcceptance/" + state + "/";
                if (!fs::exists(directory_deltaK)){
                    fs::create_directories(directory_deltaK);
                }
                std::map<int, double> trimmed_deltaK;

                for (int deltaK=0; deltaK<networkSize; ++deltaK){
                    if (sampled_deltaUpperBound_deltaAcceptance[state][deltaK] && deltaUpperBound_deltaAcceptance[state][deltaK]){
                        trimmed_deltaK[deltaK] = deltaUpperBound_deltaAcceptance[state][deltaK] / sampled_deltaUpperBound_deltaAcceptance[state][deltaK];
                    }
                }
                CSV::write(directory_deltaK + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_deltaK);
            }
        }
    }//* End of function mBFW::generate::save
}//* End of namespace mBFW::generate