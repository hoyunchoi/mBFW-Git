#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <tuple>
#include <vector>

#include "../library/CSV.hpp"
#include "../library/linearAlgebra.hpp"
#include "../library/pcg_random.hpp"
#include "NZ_Network.hpp"
#include "fileName.hpp"
#include "parameters.hpp"

namespace mBFW::generate {
const std::string baseDirectory = "../data/mBFW/";

//*--------------------------------------------Declaration of variables---------------------------------------------------------
//* Declaration of variables used at mBFW::generate namespace
int networkSize;
double acceptanceThreshold;
unsigned ensembleSize;
int coreNum;
int randomEngineSeed;
constexpr double precision = 1e4;
int maxTime;
int maxTrialTime;
const std::vector<std::string> states = {"before", "during", "after"};
double t_a, m_a, t_c, m_c;

//* Declaration of Random Engine
pcg32 randomEngine;
std::uniform_int_distribution<int> nodeDistribution;

//* Declaration of observables
//! orderParameter[time]: Average value of order parameter at specific time
// std::vector<double> orderParameter;
// std::vector<double> orderParameter_trial;

//! secondMoment[time]: Average value of second moment at specific time
// std::vector<double> secondMoment;
// std::vector<double> secondMoment_trial;

//! meanClusterSize[time]: Average value of normalized mean cluster size at specific time
// std::vector<double> meanClusterSize;
// std::vector<double> meanClusterSize_trial;

//! interEventTime[time]: Average value of inter event time at specific time
// std::vector<double> interEventTime;
// std::vector<int> sampled_interEventTime;

//! clusterSizeDist[op]: Average distribution of cluster size when order parameter passes op
//! clusterSizeDist_exact[op]: Average distribution of cluster size when order parameter is exactly op
//! clusterSizeDist_time[time]: Average distribution of cluster size when time at specific time
// std::set<double> orderParameter_clusterSizeDist;
// std::set<double> time_clusterSizeDist;
// std::map<double, std::vector<long long>> clusterSizeDist;
// std::map<double, std::vector<long long>> clusterSizeDist_exact;
// std::map<double, std::vector<long long>> clusterSizeDist_time;

//! orderParameterDist[time] : Distribution of order parameter at specific time
// std::set<double> time_orderParameterDist;
// std::map<double, std::vector<int>> orderParameterDist;

//! interEventTimeDist_op(time)["before"] : Average distribution of inter event time before discontinuous jump by order paramter (time)
//! interEventTimeDist_op(time)["during"] : Average distribution of inter event time during discontinuous jump by order paramter (time)
//! interEventTimeDist_op(time)["after"] : Average distribution of inter event time after discontinuous jump by order paramter (time)
// std::map<std::string, std::vector<int>> interEventTimeDist_op;
// std::map<std::string, std::vector<int>> interEventTimeDist_time;
// std::vector<int> interEventTimeDist_tot;

//! ageDist_op(time)["state"] : Average distribution of age before(during, after) discontinuous jump by order paramter (time)
// std::map<std::string, std::vector<long long>> ageDist_op;
// std::map<std::string, std::vector<long long>> ageDist_time;

//! deltaUpperBoundDist_op(time)["state"] : Average distribution of deltaK before(during, after) discontinuous jump by order paramter (time)
// std::map<std::string, std::vector<int>> deltaUpperBoundDist_op;
// std::map<std::string, std::vector<int>> deltaUpperBoundDist_time;
// std::vector<int> deltaUpperBoundDist_tot;

//! upperBound_deltaAcceptance[state][K] : Average delta acceptance at specific upper bound K before(during, after) discontinuous jump by order parameter
// std::map<std::string, std::vector<double>> upperBound_deltaAcceptance;
// std::map<std::string, std::vector<int>> sampled_upperBound_deltaAcceptance;

//! deltaUpperBound_deltaAcceptance[state][deltaK] : Average delta acceptance at specific delta upper bound deltaK before(during, after) discontinuous jump by order parameter
// std::map<std::string, std::vector<double>> deltaUpperBound_deltaAcceptance;
// std::map<std::string, std::vector<int>> sampled_deltaUpperBound_deltaAcceptance;

//! interEventTime_deltaUpperBound[state][iet] : Average delta upper bound at specific inter event time. log binned
// std::map<std::string, std::vector<long long>> interEventTime_deltaUpperBound;
// std::map<std::string, std::vector<int>> sampled_interEventTime_deltaUpperBound;

//! Sampled_X_interEventTime
// std::map<std::pair<int, int>, int> sampled_deltaUpperBound_interEventTime;
// std::map<std::pair<int, int>, int> sampled_upperBound_interEventTime;
// std::map<std::pair<int, int>, int> sampled_time_interEventTime;

//! noRestriction
std::vector<unsigned> noRestriction;

//! Dynacmis of network
// std::vector<std::vector<long long>> dynamics;
// std::vector<std::vector<long long>> periodDynamics;

//* Set parameters for single run
void setParameters();

//* Generate mBFW model
void run();

//* Save the model
void save();
}  // namespace mBFW::generate

void mBFW::generate::setParameters(){
    // orderParameter_clusterSizeDist = mBFW::parameters::set_orderParameter_clusterSizeDist(t_networkSize, t_acceptanceThreshold);
    // time_clusterSizeDist = mBFW::parameters::set_time_clusterSizeDist(t_networkSize, t_acceptanceThreshold);
    // time_orderParameterDist = mBFW::parameters::set_time_orderParameterDist(t_networkSize, t_acceptanceThreshold);
    std::tie(t_a, m_a, t_c, m_c) = mBFW::parameters::set_points(networkSize, acceptanceThreshold);
    t_a *= networkSize;
    t_c *= networkSize;
    maxTime = networkSize;
    maxTrialTime = std::floor(maxTime / acceptanceThreshold);

    //! Initialize Random Engine
    randomEngineSeed == -1 ? randomEngine.seed((std::random_device())()) : randomEngine.seed(randomEngineSeed);
    nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, networkSize - 1));

    //! Resize observables
    // orderParameter.assign(maxTime, 0.0);
    // secondMoment.assign(maxTime, 0.0);
    // meanClusterSize.assign(maxTime, 0.0);
    // orderParameter_trial.assign(maxTrialTime, 0.0);
    // secondMoment_trial.assign(maxTrialTime, 0.0);
    // meanClusterSize_trial.assign(maxTrialTime, 0.0);
    // interEventTime.assign(maxTime, 0.0);
    // sampled_interEventTime.assign(maxTime, 0);
    // interEventTimeDist_tot.assign(maxTime, 0);
    // deltaUpperBoundDist_tot.assign(maxTime, 0);
    noRestriction.assign(maxTime, 0);
    // for (const double& op : orderParameter_clusterSizeDist){
    //     clusterSizeDist[op].assign(t_networkSize, 0);
    // clusterSizeDist_exact[op].assign(t_networkSize, 0);
    // }
    // for (const double& t : time_clusterSizeDist){
    //     clusterSizeDist_time[t].assign(t_networkSize, 0);
    // }
    // for (const double& t : time_orderParameterDist){
    //     orderParameterDist[t].assign(t_networkSize, 0);
    // }
    // for (const std::string& state : states){
    // interEventTimeDist_op[state].assign(t_networkSize, 0);
    // interEventTimeDist_time[state].assign(t_networkSize, 0);
    // ageDist_op[state].assign(t_networkSize, 0);
    // ageDist_time[state].assign(t_networkSize, 0);
    // deltaUpperBoundDist_op[state].assign(t_networkSize, 0);
    // deltaUpperBoundDist_time[state].assign(t_networkSize, 0);
    // upperBound_deltaAcceptance[state].assign(t_networkSize, 0.0);
    // sampled_upperBound_deltaAcceptance[state].assign(t_networkSize, 0);
    // deltaUpperBound_deltaAcceptance[state].assign(t_networkSize, 0.0);
    // sampled_deltaUpperBound_deltaAcceptance[state].assign(t_networkSize, 0);
    // interEventTime_deltaUpperBound[state].assign(t_networkSize, 0);
    // sampled_interEventTime_deltaUpperBound[state].assign(t_networkSize, 0);
    // }
    // dynamics.reserve(maxTrialTime);
    // periodDynamics.reserve(maxTrialTime);
}//* End of function mBFW::generate::setParameters

void mBFW::generate::run() {
    for (unsigned ensemble = 0; ensemble < ensembleSize; ++ensemble) {
        //* Default values for one ensemble
        NZ_Network model(networkSize);
        int root1, root2;
        std::string currentState_time = "before";
        std::string currentState_op = "before";
        int size1, size2;
        int time = 0;
        int trialTime = 0;
        int periodTime = 0;
        int periodTrialTime = 0;
        int upperBound = 2;
        int eventTime = 0;
        double maxDeltaAcceptance = 0.0;
        bool choosingNewNode = true;

        //* Set initial conditions
        // orderParameter[0] += 1.0/networkSize;
        // secondMoment[0] += std::pow(1.0/networkSize, 2.0);
        // meanClusterSize[0] += 1.0;
        // orderParameter_trial[0] += 1.0/networkSize;
        // secondMoment_trial[0] += std::pow(1.0/networkSize, 2.0);
        // meanClusterSize_trial[0] += 1.0;
        // std::set<double> findingclusterSizeDist = orderParameter_clusterSizeDist;
        // std::set<double> newFindingClusterSizeDist = findingClusterSizeDist;

        //* Do mBFW algorithm until all clusters merge to one
        while (model.getMaximumClusterSize() < networkSize) {
            //* Find new nodes
            if (choosingNewNode) {
                //* Randomly choose new nodes
                do {
                    root1 = model.getRoot(nodeDistribution(randomEngine));
                    root2 = model.getRoot(nodeDistribution(randomEngine));
                } while (root1 == root2);
                //* choose two clusters of each node
                size1 = model.getClusterSize(root1);
                size2 = model.getClusterSize(root2);
            }

            //* Merge two clusters, update time
            if (size1 + size2 <= upperBound) {
                model.merge(root1, root2);
                ++time;
                ++trialTime;
                ++periodTime;
                ++periodTrialTime;
                choosingNewNode = true;
                const int currentMaximumClusterSize = model.getMaximumClusterSize();
                const double currentOrderParameter = (double)currentMaximumClusterSize / networkSize;

                //* Check the state distinguished by time
                if (time < t_a) {
                    currentState_time = "before";
                } else if (time < t_c) {
                    currentState_time = "during";
                } else {
                    currentState_time = "after";
                }

                //* Update max delta acceptance
                maxDeltaAcceptance = std::max((double)time / trialTime - acceptanceThreshold, maxDeltaAcceptance);

                //! Dynamics
                {
                    // dynamics.emplace_back(std::vector<long long> {trialTime, time, upperBound});
                    // periodDynamics.emplace_back(std::vector<long long> {periodTrialTime, periodTime, upperBound});
                }

                //! Order Parameter
                {
                    // orderParameter[time] += currentOrderParameter;
                    // orderParameter_trial[trialTime] += currentOrderParameter;
                }

                //! Second Moment
                {
                    // secondMoment[time] += std::pow(currentOrderParameter, 2.0);
                    // secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);
                }

                //! Mean Cluster Size
                {
                    // meanClusterSize[time] += model.getMeanClusterSize();
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
                    // const std::vector<std::pair<int, int>> changedAge = model.getChangedAge();
                    // for (const auto& age : changedAge){
                    //     ageDist_time[currentState_time][age.first] += age.second;
                    //     ageDist_op[currentState_op][age.first] += age.second;
                    // }
                }

                //* Order Parameter of network is changed <=> upper bound is changed
                if (model.getDeltaMaximumClusterSize() && currentMaximumClusterSize > 2) {
                    const int deltaMaximumClusterSize = model.getDeltaMaximumClusterSize();
                    //* Check the state distinguished by order parameter
                    if (currentOrderParameter < m_a) {
                        currentState_op = "before";
                    } else if (currentOrderParameter < m_c) {
                        currentState_op = "during";
                    } else {
                        currentState_op = "after";
                    }

                    //! Cluster Size Distribution
                    {
                        // findingClusterSizeDist = newFindingClusterSizeDist;
                        // for (const double& op : findingClusterSizeDist){
                        //     if (op < currentOrderParameter){
                        //         const std::map<int, int> sortedCluster = model.getSortedCluster();
                        //         for (auto it2 = sortedCluster.begin(); it2 != sortedCluster.end(); ++it2){
                        //             clusterSizeDist[op][it2->first] += it2->second;
                        //         }
                        //         newFindingClusterSizeDist.erase(op);
                        //     }
                        // }
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
                        // interEventTime[time] += (double)time-eventTime;
                        // ++sampled_interEventTime[time];
                    }

                    //! Inter Event Time Distribution time(op)
                    {
                        // ++interEventTimeDist_time[currentState_time][time-eventTime];
                        // ++interEventTimeDist_op[currentState_op][time-eventTime];
                        // ++interEventTimeDist_tot[time-eventTime];

                    }

                    //! Inter Event Time vs Delta Upper Bound
                    {
                        // interEventTime_deltaUpperBound[currentState_op][time - eventTime] += deltaMaximumClusterSize;
                        // ++sampled_interEventTime_deltaUpperBound[currentState_op][time-eventTime];
                    }

                    //! Delta Upper Bound Distribution time(op)
                    {
                        // ++deltaUpperBoundDist_time[currentState_time][deltaMaximumClusterSize];
                        // ++deltaUpperBoundDist_op[currentState_op][deltaMaximumClusterSize];
                        // ++deltaUpperBoundDist_tot[deltaMaximumClusterSize];
                    }

                    //! (Delta)Upper Bound_Delta Acceptance
                    {
                        // upperBound_deltaAcceptance[currentState_op][upperBound] += maxDeltaAcceptance;
                        // ++sampled_upperBound_deltaAcceptance[currentState_op][upperBound];
                        // deltaUpperBound_deltaAcceptance[currentState_op][deltaMaximumClusterSize] += maxDeltaAcceptance;
                        // ++sampled_deltaUpperBound_deltaAcceptance[currentState_op][upperBound];
                    }

                    //! Sampled_X_interEventTime
                    {
                        // ++sampled_deltaUpperBound_interEventTime[std::make_pair(deltaMaximumClusterSize, time - eventTime)];
                        // ++sampled_upperBound_interEventTime[std::make_pair(upperBound, time - eventTime)];
                        // ++sampled_time_interEventTime[std::make_pair(time, time - eventTime)];
                    }

                    //* Initialize variable for new period
                    eventTime = time;
                    maxDeltaAcceptance = 0;
                }  //* End of order parameter update, End of k-period
            }

            //* Upper Bound change
            else if ((double)time / trialTime <= acceptanceThreshold) {
                upperBound = size1 + size2;
                choosingNewNode = false;

                //! noRestriction
                {
                    if (periodTrialTime == 2){
                        ++noRestriction[time];
                    }
                }

                periodTime = 0;
                periodTrialTime = 0;
            }  //* End of Upper Bound change

            //* Chosen link rejected
            else {
                ++trialTime;
                ++periodTrialTime;
                choosingNewNode = true;
                // const double currentOrderParameter = (double)model.getMaximumClusterSize()/networkSize;

                //! Dynamics
                {
                    // dynamics.emplace_back(std::vector<long long> {trialTime, time, upperBound});
                    // periodDynamics.emplace_back(std::vector<long long> {periodTrialTime, periodTime, upperBound});
                }

                //! Trial Order Parameter
                // orderParameter_trial[trialTime] += currentOrderParameter;

                //! Trial Second Moment
                // secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);

                //! Trial Mean Cluster Size
                // meanClusterSize_trial[trialTime] += model.getMeanClusterSize();
            }  //* End of rejecting link
            //* End of one step
        }  //* End of network growing (one ensemble)
    }      //* End of every ensembles
}  //* End of function mBFW::generate::run


void mBFW::generate::save() {
    using namespace linearAlgebra;
    namespace fs = std::filesystem;

    //! Order Parameter
    {
        // orderParameter /= ensembleSize;
        // const std::string directory = baseDirectory + "orderParameter/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameter);
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
        // secondMoment /= ensembleSize;
        // std::vector<double> orderParameterVariance(maxTime, 0.0);
        // for (int i=0; i<maxTime; ++i){
        //     const double var = secondMoment[i] - pow(orderParameter[i], 2.0);
        //     var < 0 ? orderParameterVariance[i] = 0.0 :orderParameterVariance[i] = std::sqrt(var) * networkSize;
        // }
        // const std::string directory = baseDirectory + "orderParameterVariance/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), orderParameterVariance);
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
        // meanClusterSize /= ensembleSize;
        // const std::string directory = baseDirectory + "meanClusterSize/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), meanClusterSize);
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
        // const std::string directory = baseDirectory + "clusterSizeDist/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // for (const double& op : orderParameter_clusterSizeDist){
        //     std::map<int, double> trimmed;
        //     const long long tot = std::accumulate(clusterSizeDist[op].begin(),clusterSizeDist[op].end(), 0);
        //     for (int cs=0; cs<networkSize; ++cs){
        //         if (clusterSizeDist[op][cs]){
        //             trimmed[cs] = (double)clusterSizeDist[op][cs] / tot;
        //         }
        //     }
        //     CSV::write(directory + fileName::NGEOP(networkSize, acceptanceThreshold, ensembleSize, op, coreNum), trimmed);
        // }
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
        // const std::string directory = baseDirectory + "interEventTime/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // for (int t=0; t<maxTime; ++t){
        //     if (sampled_interEventTime[t]){
        //         interEventTime[t] /= (double)sampled_interEventTime[t];
        //     }
        // }
        // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), interEventTime);
    }

    //! Inter Event Time Distribution time(op)
    {
        // for (const std::string& state : states){
        //     const std::string directory_time = baseDirectory + "interEventTimeDist_time/" + state + "/";
        //     if (!fs::exists(directory_time)){
        //         fs::create_directories(directory_time);
        //     }
        //     std::map<int, double> trimmed_time;
        //     const double tot_time = std::accumulate(interEventTimeDist_time[state].begin(), interEventTimeDist_time[state].end(), 0.0);
        //     for (int iet=0; iet<networkSize; ++iet){
        //         if (interEventTimeDist_time[state][iet]){
        //             trimmed_time[iet] = interEventTimeDist_time[state][iet] / tot_time;
        //         }
        //     }
        //     CSV::write(directory_time + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_time);

        //     const std::string directory_op = baseDirectory + "interEventTimeDist_op/" + state + "/";
        //     if (!fs::exists(directory_op)){
        //         fs::create_directories(directory_op);
        //     }
        //     std::map<int, double> trimmed_op;
        //     const double tot_op = std::accumulate(interEventTimeDist_op[state].begin(), interEventTimeDist_op[state].end(), 0.0);
        //     for (int iet=0; iet<networkSize; ++iet){
        //         if (interEventTimeDist_op[state][iet]){
        //             trimmed_op[iet] = interEventTimeDist_op[state][iet] / tot_op;
        //         }
        //     }
        //     CSV::write(directory_op + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_op);
        // }
        // const std::string directory = baseDirectory + "interEventTimeDist_tot/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // std::map<int, double> trimmed;
        // const double tot = std::accumulate(interEventTimeDist_tot.begin(), interEventTimeDist_tot.end(), 0.0);
        // for (int iet=0; iet<networkSize; ++iet){
        //     if (interEventTimeDist_tot[iet]){
        //         trimmed[iet] = interEventTimeDist_tot[iet] / tot;
        //     }
        // }
        // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed);
    }

    //! Inter Event Time vs Delta Upper Bound
    {
        // for (const auto& state : states){
        //     const std::string directory = baseDirectory + "interEventTime_deltaUpperBound/" + state + "/";
        //     if (!fs::exists(directory)){
        //         fs::create_directories(directory);
        //     }
        //     std::map<int, double> trimmed;
        //     for (int iet=0; iet<networkSize; ++iet){
        //         if (sampled_interEventTime_deltaUpperBound[state][iet]){
        //             trimmed[iet] = (double)interEventTime_deltaUpperBound[state][iet] / sampled_interEventTime_deltaUpperBound[state][iet];
        //         }
        //     }

        //     CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed);
        // }
    }

    //! Age Distribution time(op)
    {
        // for (const auto& state : states){
        //     const std::string directory_time = baseDirectory + "ageDist_time/" + state + "/";
        //     if (!fs::exists(directory_time)){
        //         fs::create_directories(directory_time);
        //     }
        //     std::vector<double> temp_time(networkSize, 0.0);
        //     const double tot_time = std::accumulate(ageDist_time[state].begin(), ageDist_time[state].end(), 0.0);
        //     for (int age=0; age<networkSize; ++age){
        //         temp_time[age] = ageDist_time[state][age] / tot_time;
        //     }
        //     print(temp_time);
        //     CSV::write(directory_time + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), temp_time);

        //     const std::string directory_op = baseDirectory + "ageDist_op/" + state + "/";
        //     if (!fs::exists(directory_op)){
        //         fs::create_directories(directory_op);
        //     }
        //     std::vector<double> temp_op(networkSize, 0.0);
        //     const double tot_op = std::accumulate(ageDist_op[state].begin(), ageDist_op[state].end(), 0.0);
        //     for (int age=0; age<networkSize; ++age){
        //         temp_op[age] = ageDist_time[state][age] / tot_time;
        //     }
        //     CSV::write(directory_op + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), temp_op);
        // }
    }

    //! Delta Upper Bound Distribution_time(op)
    {
        // for (const auto& state : states){
        //     const std::string directory_time = baseDirectory + "deltaUpperBoundDist_time/" + state + "/";
        //     if (!fs::exists(directory_time)){
        //         fs::create_directories(directory_time);
        //     }
        //     std::map<int, double> trimmed_time;
        //     const double tot_time = std::accumulate(deltaUpperBoundDist_time[state].begin(), deltaUpperBoundDist_time[state].end(), 0.0);
        //     for (int deltaK = 0; deltaK < networkSize; ++deltaK){
        //         if (deltaUpperBoundDist_time[state][deltaK]){
        //             trimmed_time[deltaK] = (double)deltaUpperBoundDist_time[state][deltaK] / tot_time;
        //         }
        //     }
        //     CSV::write(directory_time + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_time);

        //     const std::string directory_op = baseDirectory + "deltaUpperBoundDist_op/" + state + "/";
        //     if (!fs::exists(directory_op)){
        //         fs::create_directories(directory_op);
        //     }
        //     std::map<int, double> trimmed_op;
        //     const double tot_op = std::accumulate(deltaUpperBoundDist_op[state].begin(), deltaUpperBoundDist_op[state].end(), 0.0);
        //     for (int deltaK=0; deltaK<networkSize; ++deltaK){
        //         if (deltaUpperBoundDist_op[state][deltaK]){
        //             trimmed_op[deltaK] = (double)deltaUpperBoundDist_op[state][deltaK] / tot_op;
        //         }
        //     }
        //     CSV::write(directory_op + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_op);
        // }

        // const std::string directory = baseDirectory + "deltaUpperBoundDist_tot/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // std::map<int, double> trimmed;
        // const double tot = std::accumulate(deltaUpperBoundDist_tot.begin(), deltaUpperBoundDist_tot.end(), 0.0);
        // for (int deltaK = 0; deltaK < networkSize; ++deltaK){
        //     if (deltaUpperBoundDist_tot[deltaK]){
        //         trimmed[deltaK] = (double)deltaUpperBoundDist_tot[deltaK] / tot;
        //     }
        // }
        // CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed);
    }

    //! (Delta) Upper Bound vs Delta Acceptance
    {
        // for (const auto& state : states){
        //     const std::string directory_k = baseDirectory + "upperBound_deltaAcceptance/" + state + "/";
        //     if (!fs::exists(directory_k)){
        //         fs::create_directories(directory_k);
        //     }
        //     std::map<int, double> trimmed_k;
        //     for (int k=0; k<networkSize; ++k){
        //         if (sampled_upperBound_deltaAcceptance[state][k] && upperBound_deltaAcceptance[state][k]){
        //             trimmed_k[k] = upperBound_deltaAcceptance[state][k] / sampled_upperBound_deltaAcceptance[state][k];
        //         }
        //     }
        //     CSV::write(directory_k + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_k);

        //     const std::string directory_deltaK = baseDirectory + "deltaUpperBound_deltaAcceptance/" + state + "/";
        //     if (!fs::exists(directory_deltaK)){
        //         fs::create_directories(directory_deltaK);
        //     }
        //     std::map<int, double> trimmed_deltaK;

        //     for (int deltaK=0; deltaK<networkSize; ++deltaK){
        //         if (sampled_deltaUpperBound_deltaAcceptance[state][deltaK] && deltaUpperBound_deltaAcceptance[state][deltaK]){
        //             trimmed_deltaK[deltaK] = deltaUpperBound_deltaAcceptance[state][deltaK] / sampled_deltaUpperBound_deltaAcceptance[state][deltaK];
        //         }
        //     }
        //     CSV::write(directory_deltaK + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_deltaK);
        // }
    }

    //! Sampled_X_interEventTime
    {
        // const std::string directory_dk = baseDirectory + "sampled_deltaUpperBound_interEventTime/";
        // if (!fs::exists(directory_dk)) {
        //     fs::create_directories(directory_dk);
        // }
        // std::map<std::pair<int, int>, double> trimmed_dk(sampled_deltaUpperBound_interEventTime.begin(), sampled_deltaUpperBound_interEventTime.end());
        // const double tot_dk = accumulate(trimmed_dk);
        // trimmed_dk /= tot_dk;
        // CSV::write(directory_dk + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_dk);

        // const std::string directory_k = baseDirectory + "sampled_upperBound_interEventTime/";
        // if (!fs::exists(directory_k)) {
        //     fs::create_directories(directory_k);
        // }
        // std::map<std::pair<int, int>, double> trimmed_k(sampled_upperBound_interEventTime.begin(), sampled_upperBound_interEventTime.end());
        // const double tot_k = accumulate(trimmed_k);
        // trimmed_k /= tot_k;
        // CSV::write(directory_k + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_k);

        // const std::string directory_t = baseDirectory + "sampled_time_interEventTime/";
        // if (!fs::exists(directory_t)) {
        //     fs::create_directories(directory_t);
        // }
        // std::map<std::pair<int, int>, double> trimmed_t(sampled_time_interEventTime.begin(), sampled_time_interEventTime.end());
        // const double tot_t = accumulate(trimmed_t);
        // trimmed_t /= tot_t;
        // CSV::write(directory_t + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed_t);
    }

    //! No Restriction
    {
        const std::string directory = baseDirectory + "noRestriction/";
        if (!fs::exists(directory)){
            fs::create_directories(directory);
        }
        std::map<int, double> trimmed;
        for (int t=0; t<maxTime; ++t){
            if (noRestriction[t]){
                trimmed[t] = noRestriction[t] / (double)ensembleSize;
            }
        }
        CSV::write(directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum), trimmed);
    }

    //! Dynamics
    {
        // const std::string directory = baseDirectory + "dynamics/";
        // if (!fs::exists(directory)){
        //     fs::create_directories(directory);
        // }
        // CSV::write(directory + fileName::base(networkSize, acceptanceThreshold) + "-" + std::to_string(randomEngineSeed) + ".txt", dynamics);
        // const std::string periodDirectory = baseDirectory + "periodDynamics/";
        // if (!fs::exists(periodDirectory)){
        //     fs::create_directories(periodDirectory);
        // }
        // CSV::write(periodDirectory + fileName::base(networkSize, acceptanceThreshold) + "-" + std::to_string(randomEngineSeed) + ".txt", periodDynamics);
    }
}  //* End of function mBFW::generate::save