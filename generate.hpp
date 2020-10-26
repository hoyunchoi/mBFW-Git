#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>

#include "../library-Git/linearAlgebra.hpp"
#include "../library-Git/Networks.hpp"
#include "../library-Git/CSV.hpp"
#include "../library-Git/pcg_random.hpp"

#include "parameters.hpp"
#include "mBFW.hpp"

namespace mBFW::generate{
    //*--------------------------------------------Declaration of variables---------------------------------------------------------
    //! Declaration of variables only used at mBFW::generate namespace
    int randomEngineSeed;
    double precision;
    double degenerated;

    //* variables for log binning
    int logBinNum;
    std::vector<double> minBin;
    std::vector<double> valueBin;

    //! Random Engine
    pcg32 randomEngine;
    std::uniform_int_distribution<int> nodeDistribution;
    int a;

    //! calculation functions
    std::function<void(const int&, const double&)> do_orderParameter;
    std::function<void(const int&, const double&)> do_meanClusterSize;
    std::function<void(const int&, const double&)> do_secondGiant;
    std::function<void(const std::string&, const int&, const int&)> do_interEventTime;
    std::function<void(const std::string&, const int&, const double&)> do_deltaAcceptance;
    std::function<void(const int&, const int&)> do_orderParameterDistribution;
    std::function<void(const double&, NZ_Network&)> do_clusterSizeDistribution;
    std::function<void(const double&, const NZ_Network&)> do_ageDistribution;
    std::function<void(const std::string&, const int&)> do_interEventTimeDistribution;
    std::function<void(const std::string&, const int&)> do_deltaUpperBoundDistribution;
    std::function<void(const std::string&, const double&)> do_deltaAcceptanceDistribution;
    std::function<void(const int&, const double&)> do_interEventTime_DeltaAcceptance;
    std::function<void(const int&, const double&)> do_upperBound_DeltaAcceptance;
    std::function<void(const int&, const double&)> do_deltaUpperBound_DeltaAcceptance;
    std::function<void(const double&, const int&, const int&, const int&, const int&, const int&, const int&)> do_dynamics;

    //! Declaration of observables
    //* X[time] = value of observable X at time
    std::vector<double> orderParameter;
    std::vector<double> secondGiant;
    std::vector<double> meanClusterSize;
    std::map<std::string, std::vector<double>> interEventTime;
    std::map<std::string, std::vector<int>> sampledInterEventTime;
    std::map<std::string, std::vector<double>> deltaAcceptance;
    std::map<std::string, std::vector<int>> sampledDeltaAcceptance;

    //* orderParameterDistribution[time] : distribution of rounded order parameter at (rounded time)=(time)
    //* orderParameterDistribution[time][mcs] : number of cluster size of "maximum cluster size"
    std::map<double, std::vector<int>> orderParameterDistribution;

    //* clusterSizeDistribution[orderParameter] : cluster size distribution at (rounded order parameter)=(orderParameter)
    //* clusterSizeDistribution[orderParameter][cs] : number of cluster of size "cs"
    std::map<double, std::vector<long long>> clusterSizeDistribution;
    std::map<double, long long> sampledClusterSizeDistribution;

    //* interEventTimeDistribution["before"] : inter event time distribution before jump (order parameter < m_a)
    //* interEventTimeDistribution["during"] : inter event time distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<int>> interEventTimeDistribution;

    //* detaMDistribution["before"] : jump size of maximum cluster size distribution before jump
    //* detaMDistribution["during"] : jump size of maximum cluster size distribution during jump
    std::map<std::string, std::vector<int>> deltaUpperBoundDistribution;

    //* ageDistribution["before"] : age distribution before jump (order parameter < m_a)
    //* ageDistribution["during"] : age distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<long long>> ageDistribution;

    //* deltaAcceptanceDistribution["before"] : max acceptance distribution befeore jump (order paramter < m_a)
    //* deltaAcceptanceDistribution["during"] : max acceptance distribution during jump (m_a < order paramter < m_c)
    std::map<std::string, std::vector<int>> deltaAcceptanceDistribution;

    //* X_maxAcceptance[x] : average max acceptance at X=x before jump
    //* sampledX_maxAcceptance[x] : number of samples for each X=x
    std::vector<double> interEventTime_DeltaAcceptance;
    std::vector<int> sampledInterEventTime_DeltaAcceptance;
    std::vector<double> upperBound_DeltaAcceptance;
    std::vector<int> sampledUpperBound_DeltaAcceptance;
    std::vector<double> deltaUpperBound_DeltaAcceptance;
    std::vector<int> sampledDeltaUpperBound_DeltaAcceptance;

    //* Dynamics
    std::map<std::string, std::vector<std::vector<int>>> dynamics;

    //*------------------------------------------- functions for calculating observables ------------------------------------------------------
    void calculate_orderParameter(const int& t_time, const double& t_exactOrderParameter){
        orderParameter[t_time] += t_exactOrderParameter;
    }
    void calculate_meanClusterSize(const int& t_time, const double& t_meanClusterSize){
        meanClusterSize[t_time] += t_meanClusterSize;
    }
    void calculate_secondGiant(const int& t_time, const double& t_secondMaximumCluster){
        secondGiant[t_time] += t_secondMaximumCluster;
    }
    void calculate_interEventTime(const std::string& t_currentState, const int& t_time, const int& t_currentInterEventTime){
        interEventTime[t_currentState][t_time] += t_currentInterEventTime;
        ++sampledInterEventTime[t_currentState][t_time];
    }
    void calculate_deltaAcceptance(const std::string& t_currentState, const int& t_time, const double& t_currentDeltaAcceptance){
        deltaAcceptance[t_currentState][t_time] += t_currentDeltaAcceptance;
        ++sampledDeltaAcceptance[t_currentState][t_time];
    }
    void calculate_orderParameterDistribution(const int& t_time, const int& t_maximumClusterSize){
        const double roundedTime = round(t_time/degenerated)/precision;
        auto it = std::find(time_orderParameterDistribution.begin(), time_orderParameterDistribution.end(), roundedTime);
        if (it != time_orderParameterDistribution.end()){
            ++orderParameterDistribution[*it][t_maximumClusterSize];
        }
    }
    void calculate_clusterSizeDistribution(const double& t_exactOrderParameter, NZ_Network& t_model){
        const double roundedOrderParameter = round(t_exactOrderParameter*precision)/precision;
        auto it = std::find(orderParameter_clusterSizeDistribution.begin(), orderParameter_clusterSizeDistribution.end(), roundedOrderParameter);
        if (it != orderParameter_clusterSizeDistribution.end()){
            ++sampledClusterSizeDistribution[*it];
            const std::map<int,int> sortedCluster = t_model.getSortedCluster();
            for (auto it2 = sortedCluster.begin(); it2!= sortedCluster.end(); ++it2){
                clusterSizeDistribution[*it][it2->first] += it2->second;
            }
        }
    }
    void calculate_ageDistribution(const double& t_exactOrderParameter, const NZ_Network& t_model){
        if (t_exactOrderParameter < m_c){
            const std::vector<std::pair<int,int>> changedAge = t_model.getChangedAge();
            std::string currentState;
            t_exactOrderParameter < m_a ? currentState = "before" : currentState = "during";
            for (const auto& age : changedAge){
                ageDistribution[currentState][age.first] += age.second;
            }
        }
    }
    void calculate_interEventTimeDistribution(const std::string& t_currentState, const int& t_currentInterEventTime){
        ++interEventTimeDistribution[t_currentState][t_currentInterEventTime];
    }
    void calculate_deltaUpperBoundDistribution(const std::string& t_currentState, const int& t_deltaMaximumClusterSize){
        ++deltaUpperBoundDistribution[t_currentState][t_deltaMaximumClusterSize];
    }
    void calculate_deltaAcceptanceDistribution(const std::string& t_currentState, const double& t_currentDeltaAcceptance){
        for (int i=0; i<logBinNum; ++i){
            if (minBin[i+1] > t_currentDeltaAcceptance && t_currentDeltaAcceptance){
                ++deltaAcceptanceDistribution[t_currentState][i];
                break;
            }
        }
    }
    void calculate_interEventTime_DeltaAcceptance(const int& t_currentInterEventTime, const double& t_currentDeltaAcceptance){
        interEventTime_DeltaAcceptance[t_currentInterEventTime] += t_currentDeltaAcceptance;
        ++sampledInterEventTime_DeltaAcceptance[t_currentInterEventTime];
    }
    void calculate_upperBound_DeltaAcceptance(const int& t_upperBound, const double& t_currentDeltaAcceptance){
        upperBound_DeltaAcceptance[t_upperBound] += t_currentDeltaAcceptance;
        ++sampledUpperBound_DeltaAcceptance[t_upperBound];}
    void calculate_deltaUpperBound_DeltaAcceptance(const int& t_deltaMaximumClusterSize, const double& t_currentDeltaAcceptance){
        deltaUpperBound_DeltaAcceptance[t_deltaMaximumClusterSize] += t_currentDeltaAcceptance;
        ++sampledDeltaUpperBound_DeltaAcceptance[t_deltaMaximumClusterSize];
    }
    void calculate_dynamics(const double& t_exactOrderParameter, const int& t_time, const int& t_trialTime, const int& t_periodTime, const int& t_periodTrialTime, const int& t_maximumClusterSize, const int& t_upperBound){
        if (t_exactOrderParameter < m_c){
            std::string currentState;
            t_exactOrderParameter < m_a ? currentState = "before" : currentState = "during";
             dynamics[currentState].emplace_back(std::vector<int>{t_time, t_trialTime, t_periodTime, t_periodTrialTime, t_maximumClusterSize, t_upperBound});
        }
    }

    //*-------------------------------------------Set Parameters for one run------------------------------------------------------
    void setParameters(const int& t_networkSize, const int& t_ensembleSize, const double& t_acceptanceThreshold, const double t_precision, const int& t_coreNum, const int& t_randomEngineSeed, const std::vector<bool> t_observables){
        //! Observables to be processed
        process_orderParameter = t_observables[0];
        process_meanClusterSize = t_observables[1];
        process_secondGiant = t_observables[2];
        process_interEventTime = t_observables[3];
        process_deltaAcceptance = t_observables[4];
        process_orderParameterDistribution = t_observables[5];
        process_clusterSizeDistribution = t_observables[6];
        process_ageDistribution = t_observables[7];
        process_interEventTimeDistribution = t_observables[8];
        process_deltaUpperBoundDistribution = t_observables[9];
        process_deltaAcceptanceDistribution = t_observables[10];
        process_interEventTime_DeltaAcceptance = t_observables[11];
        process_upperBound_DeltaAcceptance = t_observables[12];
        process_deltaUpperBound_DeltaAcceptance = t_observables[13];
        process_dynamics = t_observables[14];

        //! Input variables
        networkSize = t_networkSize;
        ensembleSize = t_ensembleSize;
        coreNum = t_coreNum;
        acceptanceThreshold = t_acceptanceThreshold;
        randomEngineSeed = t_randomEngineSeed;
        t_networkSize < t_precision ? precision = t_networkSize : precision = t_precision;

        //! Random Engine
        randomEngineSeed == -1 ? randomEngine.seed((std::random_device())()) : randomEngine.seed(randomEngineSeed);
        nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, networkSize-1));

        //! calculation function
        process_orderParameter ? do_orderParameter = calculate_orderParameter : do_orderParameter = [](const int&, const double&){};
        process_meanClusterSize ? do_meanClusterSize = calculate_meanClusterSize : do_meanClusterSize = [](const int&, const double&){};
        process_secondGiant ? do_secondGiant = calculate_secondGiant : do_secondGiant = [](const int&, const double&){};
        process_interEventTime ? do_interEventTime = calculate_interEventTime : do_interEventTime = [](const std::string&, const int&, const int&){};
        process_deltaAcceptance ? do_deltaAcceptance = calculate_deltaAcceptance : do_deltaAcceptance = [](const std::string&, const int&, const double&){};
        process_orderParameterDistribution ? do_orderParameterDistribution = calculate_orderParameterDistribution : do_orderParameterDistribution = [](const int&, const int&){};
        process_clusterSizeDistribution ? do_clusterSizeDistribution = calculate_clusterSizeDistribution : do_clusterSizeDistribution = [](const double&, NZ_Network&){};
        process_ageDistribution ? do_ageDistribution = calculate_ageDistribution : do_ageDistribution = [](const double&, const NZ_Network&){};
        process_interEventTimeDistribution ? do_interEventTimeDistribution = calculate_interEventTimeDistribution : do_interEventTimeDistribution = [](const std::string&, const int&){};
        process_deltaUpperBoundDistribution ? do_deltaUpperBoundDistribution = calculate_deltaUpperBoundDistribution : do_deltaUpperBoundDistribution = [](const std::string&, const int&){};
        process_deltaAcceptanceDistribution ? do_deltaAcceptanceDistribution = calculate_deltaAcceptanceDistribution : do_deltaAcceptanceDistribution = [](const std::string&, const double&){};
        process_interEventTime_DeltaAcceptance ? do_interEventTime_DeltaAcceptance = calculate_interEventTime_DeltaAcceptance : do_interEventTime_DeltaAcceptance = [](const int&, const double&){};
        process_upperBound_DeltaAcceptance ? do_upperBound_DeltaAcceptance = calculate_upperBound_DeltaAcceptance : do_upperBound_DeltaAcceptance = [](const int&, const double&){};
        process_deltaUpperBound_DeltaAcceptance ? do_deltaUpperBound_DeltaAcceptance = calculate_deltaUpperBound_DeltaAcceptance : do_deltaUpperBound_DeltaAcceptance = [](const int&, const double&){};
        process_dynamics ? do_dynamics = calculate_dynamics : do_dynamics = [](const double&, const int&, const int&, const int&, const int&, const int&, const int&){};

        //! Pre-defined variables
        std::tie(m_c, t_c, time_orderParameterDistribution, orderParameter_clusterSizeDistribution) = mBFW::parameters::pre_defined(networkSize, acceptanceThreshold);
        degenerated=t_networkSize/t_precision;

        //! Output variables (Observables)
        //* time-X
        if (process_orderParameter){orderParameter.assign(t_networkSize, 0.0);}
        if (process_meanClusterSize){meanClusterSize.assign(t_networkSize, 0.0);}
        if (process_secondGiant){secondGiant.assign(t_networkSize, 0.0);}
        for (auto state : states){
            if (process_interEventTime){
                interEventTime[state].assign(t_networkSize, 0.0);
                sampledInterEventTime[state].assign(t_networkSize, 0);
            }
            if (process_deltaAcceptance){
                deltaAcceptance[state].assign(t_networkSize, 0.0);
                sampledDeltaAcceptance[state].assign(t_networkSize, 0);
            }
        }

        //* Distributions
        if (process_orderParameterDistribution){
            for (const double& t : time_orderParameterDistribution){
                orderParameterDistribution[t].assign(t_networkSize, 0);
            }
        }
        if (process_clusterSizeDistribution){
            for (const double& m : orderParameter_clusterSizeDistribution){
                clusterSizeDistribution[m].assign(t_networkSize, 0);
            }
        }
        for (auto state : states){
            if (process_interEventTimeDistribution){interEventTimeDistribution[state].assign(t_networkSize, 0);}
            if (process_deltaUpperBoundDistribution){deltaUpperBoundDistribution[state].assign(t_networkSize, 0);}
            if (process_ageDistribution){ageDistribution[state].assign(t_networkSize, 0);}
            if (process_deltaAcceptanceDistribution){
                const std::vector<double> exponent = linearAlgebra::arange(-8.0, 0.0, 0.01);
                minBin = linearAlgebra::elementPow(10.0, exponent);
                logBinNum = exponent.size()-1;
                valueBin.assign(logBinNum, 0.0);
                for (int i=0; i<logBinNum; ++i){
                    valueBin[i] = sqrt(minBin[i] * minBin[i+1]);
                }
                deltaAcceptanceDistribution[state].assign(logBinNum, 0);
            }
        }

        //* X-DeltaAcceptance
        if (process_interEventTime_DeltaAcceptance){
            interEventTime_DeltaAcceptance.assign(t_networkSize, 0.0);
            sampledInterEventTime_DeltaAcceptance.assign(t_networkSize, 0);
        }
        if (process_upperBound_DeltaAcceptance){
            upperBound_DeltaAcceptance.assign(t_networkSize, 0.0);
            sampledUpperBound_DeltaAcceptance.assign(t_networkSize, 0);
        }
        if (process_deltaUpperBound_DeltaAcceptance){
            deltaUpperBound_DeltaAcceptance.assign(t_networkSize, 0.0);
            sampledDeltaUpperBound_DeltaAcceptance.assign(t_networkSize, 0);
        }
        //* Dynamics
        if (process_dynamics){
            if (t_ensembleSize < 5){
                dynamics["before"].reserve(networkSize*ensembleSize);
                dynamics["during"].reserve(networkSize*ensembleSize);
            }
            else{
                std::cout<<"To many ensemble size to save dynamics.\n";
                exit(0);
            }
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
            int periodTime = 0;
            int periodTrialTime = 0;
            int peakTime = 0;
            int peakTrialTime = 0;
            int updatedTime = 0;
            double maxAcceptance = 0;
            double maxPeriodAcceptance = 0;
            bool findNewNodes = true;

            //* initial condition
            if (process_orderParameter){orderParameter[0] += 1.0/networkSize;}
            if (process_meanClusterSize){meanClusterSize[0] += 1.0;}
            if (process_secondGiant){secondGiant[0] += 1.0/networkSize;}

            //* Do rBFW algorithm until all clusters merge to one
            while (model.getMaximumClusterSize() < networkSize){
                //* Find new nodes
                if (findNewNodes){
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
                    ++periodTime;
                    ++periodTrialTime;
                    findNewNodes=true;
                    const int maximumClusterSize = model.getMaximumClusterSize();
                    const double exactOrderParameter = maximumClusterSize/(double)networkSize;

                    //* max acceptance
                    if ((double)time/trialTime > maxAcceptance){
                        peakTime = time;
                        peakTrialTime = trialTime;
                        maxAcceptance = (double)time/trialTime;
                    }

                    //! Dynamics
                    do_dynamics(exactOrderParameter, time, trialTime, periodTime, periodTrialTime, maximumClusterSize, upperBound);

                    //! Order Parameter
                    do_orderParameter(time, exactOrderParameter);

                    //! Mean cluster Size
                    do_meanClusterSize(time, model.getMeanClusterSize());

                    //! Second Giant
                    do_secondGiant(time, model.getSecondMaximumClusterSize()/(double)networkSize);

                    //! Age Distribution
                    do_ageDistribution(exactOrderParameter, model);

                    //! Order Parameter Distribution
                    do_orderParameterDistribution(time, maximumClusterSize);

                    //* End of k-period
                    if (model.getDeltaMaximumClusterSize() && model.getMaximumClusterSize()>2){
                        const int currentInterEventTime = time-updatedTime;
                        const double currentDeltaAcceptance = maxAcceptance-acceptanceThreshold;
                        std::string currentState;

                        //! Cluster Size Distribution
                        do_clusterSizeDistribution(exactOrderParameter, model);

                        //* Before and During Jump
                        if (exactOrderParameter < m_c){
                            const int deltaMaximumClusterSize = model.getDeltaMaximumClusterSize();
                            exactOrderParameter < m_a ? currentState = "before" : currentState = "during";

                            //! Inter Event Time Distribution
                            do_interEventTimeDistribution(currentState, currentInterEventTime);

                            //! Delta Upper Bound Distribution
                            do_deltaUpperBoundDistribution(currentState, deltaMaximumClusterSize);

                            //! Delta Acceptance Distribution
                            do_deltaAcceptanceDistribution(currentState, currentDeltaAcceptance);

                            //! Inter Event Time
                            do_interEventTime(currentState, time, currentInterEventTime);

                            //! Delta Acceptance
                            do_deltaAcceptance(currentState, time, currentDeltaAcceptance);

                            //! X vs Delta Acceptance
                            if (exactOrderParameter < m_a){
                                do_interEventTime_DeltaAcceptance(currentInterEventTime, currentDeltaAcceptance);
                                do_upperBound_DeltaAcceptance(upperBound, currentDeltaAcceptance);
                                do_deltaUpperBound_DeltaAcceptance(deltaMaximumClusterSize, currentDeltaAcceptance);
                            }
                        }

                        //* Initialize variable for new k-period
                        updatedTime = time;
                        periodTime = 0;
                        periodTrialTime = 0;
                        maxAcceptance=0;
                        maxPeriodAcceptance=0;
                    }//* End of updating k-period
                }
                else if ((double)time/trialTime < acceptanceThreshold){
                    ++upperBound;
                    findNewNodes=false;
                }
                else{
                    ++trialTime;
                    ++periodTrialTime;
                    findNewNodes=true;
                }//* End of one step
            }//* End of network growing (one ensemble)
        } //* End of every ensembles
    } //* End of function mBFW::generate::run

    //*-------------------------------------------Save calculated variables------------------------------------------------------
    void save(){
        using namespace linearAlgebra;

        //! Order Parameter
        if (process_orderParameter){
            const std::string fullFileName = rootDirectory + "orderParameter/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            orderParameter /= ensembleSize;
            writeCSV(fullFileName, orderParameter);
        }

        //! Mean Cluster Size
        if (process_meanClusterSize){
            const std::string fullFileName = rootDirectory + "meanClusterSize/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            meanClusterSize /= ensembleSize;
            writeCSV(fullFileName, meanClusterSize);
        }

        //! Second Giant
        if (process_secondGiant){
            const std::string fullFileName = rootDirectory + "secondGiant/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            secondGiant /= ensembleSize;
            writeCSV(fullFileName, secondGiant);
        }

        //! Inter Event Time
        if (process_interEventTime){
            for (const auto& state : states){
                const std::string fullFileName = rootDirectory + "interEventTime/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

                //* save only useful data
                std::map<double, double> trimmed;
                for (int t=0; t<networkSize; ++t){
                    if (sampledInterEventTime[state][t]){
                        trimmed[(double)t/networkSize] = interEventTime[state][t]/sampledInterEventTime[state][t];
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Delta Acceptance
        if (process_deltaAcceptance){
            for (const auto& state : states){
                const std::string fullFileName = rootDirectory + "deltaAcceptance/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

                //* save only useful data
                std::map<double, double> trimmed;
                for (int t=0; t<networkSize; ++t){
                    if (sampledDeltaAcceptance[state][t]){
                        trimmed[(double)t/networkSize] = deltaAcceptance[state][t]/sampledDeltaAcceptance[state][t];
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Order Parameter Distribution
        if (process_orderParameterDistribution){
            for (const auto& t : time_orderParameterDistribution){
                const std::string fullFileName = rootDirectory + "orderParameterDistribution/" + filename_time(networkSize, acceptanceThreshold, ensembleSize, t, coreNum);

                //* save only useful data
                std::map<double, double> trimmed;
                const double tot = std::accumulate(orderParameterDistribution[t].begin(), orderParameterDistribution[t].end(), 0);
                for (int mcs=0; mcs<networkSize; ++mcs){
                    if (orderParameterDistribution[t][mcs]){
                        trimmed[(double)mcs/networkSize] = orderParameterDistribution[t][mcs]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Cluster Size Distribution
        if (process_clusterSizeDistribution){
            for (const auto& op : orderParameter_clusterSizeDistribution){
                const std::string fullFileName = rootDirectory + "clusterSizeDistribution/" + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSize, op, coreNum);

                //* save only useful data
                std::map<int, double> trimmed;
                if (sampledClusterSizeDistribution[op] != 0){
                    for (int cs=0; cs<networkSize; ++cs){
                        if (clusterSizeDistribution[op][cs]){
                            trimmed.insert(std::make_pair(cs, (double)clusterSizeDistribution[op][cs]/sampledClusterSizeDistribution[op]));
                        }
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Age Distribution
        if (process_ageDistribution){
            for (const auto& state : states){
                const std::string fullFileName = rootDirectory + "ageDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

                //* save only useful data
                std::map<int, double> trimmed;
                const double tot = std::accumulate(ageDistribution[state].begin(), ageDistribution[state].end(), 0.0);
                for (int age=0; age<networkSize; ++age){
                    if (ageDistribution[state][age]){
                        trimmed[age] = ageDistribution[state][age]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Inter Event Time Distribution
        if (process_interEventTimeDistribution){
            for (const auto& state : states){
                const std::string fullFileName = rootDirectory + "interEventTimeDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

                //* save only useful data
                std::map<int, double> trimmed;
                const double tot = std::accumulate(interEventTimeDistribution[state].begin(), interEventTimeDistribution[state].end(), 0.0);
                for (int iet=0; iet<networkSize; ++iet){
                    if (interEventTimeDistribution[state][iet]){
                        trimmed[iet] = interEventTimeDistribution[state][iet]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Delta Upper Bound Distribution
        if (process_deltaUpperBoundDistribution){
            for (const auto& state : states){
                const std::string fullFileName = rootDirectory + "deltaUpperBoundDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

                //* save only useful data
                std::map<int, double> trimmed;
                const double tot = std::accumulate(deltaUpperBoundDistribution[state].begin(), deltaUpperBoundDistribution[state].end(), 0.0);
                for (int deltaK=0; deltaK<networkSize; ++deltaK){
                    if (deltaUpperBoundDistribution[state][deltaK]){
                        trimmed[deltaK] = deltaUpperBoundDistribution[state][deltaK]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Delta Acceptance Distribution
        if (process_deltaAcceptanceDistribution){
            for (const auto& state : states){
                const std::string fullFileName = rootDirectory + "deltaAcceptanceDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

                //* save only useful data
                std::map<double, double> trimmed;
                const int tot = std::accumulate(deltaAcceptanceDistribution[state].begin(), deltaAcceptanceDistribution[state].end(), 0);
                for (int i=0; i<logBinNum; ++i){
                    if (deltaAcceptanceDistribution[state][i]){
                        trimmed[valueBin[i]] = (double)deltaAcceptanceDistribution[state][i]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Inter Event Time vs Delta Acceptance
        if (process_interEventTime_DeltaAcceptance){
            const std::string fullFileName = rootDirectory + "interEventTime_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

            //* save only useful data
            std::map<int, double> trimmed;
            for (int iet=1; iet<networkSize; ++iet){
                if (sampledInterEventTime_DeltaAcceptance[iet] && interEventTime_DeltaAcceptance[iet]){
                    trimmed[iet] = interEventTime_DeltaAcceptance[iet]/sampledInterEventTime_DeltaAcceptance[iet];
                }
            }
            writeCSV(fullFileName, trimmed);
        }

        //! Upper Bound vs Delta Acceptance
        if (process_upperBound_DeltaAcceptance){
            const std::string fullFileName = rootDirectory + "upperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

            //* save only useful data
            std::map<int, double> trimmed;
            for (int k=2; k<networkSize; ++k){
                if (sampledUpperBound_DeltaAcceptance[k] && upperBound_DeltaAcceptance[k]){
                    trimmed[k] = upperBound_DeltaAcceptance[k]/sampledUpperBound_DeltaAcceptance[k];
                }
            }
            writeCSV(fullFileName, trimmed);
        }

        //! Delta Upper Bound vs Delta Acceptance
        if (process_deltaUpperBound_DeltaAcceptance){
            const std::string fullFileName = rootDirectory + "deltaUpperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

            //* save only useful data
            std::map<int, double> trimmed;
            for (int deltaK=1; deltaK<networkSize; ++deltaK){
                if (sampledDeltaUpperBound_DeltaAcceptance[deltaK] && deltaUpperBound_DeltaAcceptance[deltaK]){
                    trimmed[deltaK] = deltaUpperBound_DeltaAcceptance[deltaK]/sampledDeltaUpperBound_DeltaAcceptance[deltaK];
                }
            }
            writeCSV(fullFileName, trimmed);
        }

        //! Dynamics
        if (process_dynamics){
            for (auto state : states){
                const std::string fullFileName = rootDirectory + "dynamics/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, -1, randomEngineSeed);
                writeCSV(fullFileName, dynamics[state]);

            }
        }
    }
}//* End of namespace mBFW::generate