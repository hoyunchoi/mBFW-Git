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

    //*-------------------------------------------Set Parameters for one run------------------------------------------------------
    void setParameters(const int& t_networkSize, const int& t_ensembleSize, const double& t_acceptanceThreshold, const int& t_coreNum, const int& t_randomEngineSeed){
        //* Input variables
        networkSize = t_networkSize;
        ensembleSize = t_ensembleSize;
        coreNum = t_coreNum;
        acceptanceThreshold = t_acceptanceThreshold;
        randomEngineSeed = t_randomEngineSeed;

        //* Initialize Random Engine
        randomEngineSeed == -1 ? randomEngine.seed((std::random_device())()) : randomEngine.seed(randomEngineSeed);
        nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, networkSize-1));

        //* Resize observables
        const int maxTime = t_networkSize;
        const int maxTrialTime = std::ceil(t_networkSize/t_acceptanceThreshold);
        orderParameter.assign(maxTime, 0.0);
        orderParameter_trial.assign(maxTrialTime, 0.0);
        secondMoment.assign(maxTime, 0.0);
        secondMoment_trial.assign(maxTrialTime, 0.0);
        meanClusterSize.assign(maxTime, 0.0);
        meanClusterSize_trial.assign(maxTrialTime, 0.0);
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
            bool findNewNodes = true;

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
                    findNewNodes = true;
                    const double currentOrderParameter = (double)model.getMaximumClusterSize()/networkSize;

                    //! Order Parameter
                    orderParameter[time] += currentOrderParameter;
                    orderParameter_trial[trialTime] += currentOrderParameter;

                    //! Second Moment
                    secondMoment[time] += std::pow(currentOrderParameter, 2.0);
                    secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);

                    //! Mean Cluster Size
                    meanClusterSize[time] += model.getMeanClusterSize();
                    meanClusterSize_trial[trialTime] += model.getMeanClusterSize();
                }
                else if ((double)time/trialTime <= acceptanceThreshold){
                    upperBound = size1+size2;
                    findNewNodes=false;
                }
                else{
                    ++trialTime;
                    findNewNodes=true;
                    const double currentOrderParameter = (double)model.getMaximumClusterSize()/networkSize;

                    //! Trial Order Parameter
                    orderParameter_trial[trialTime] += currentOrderParameter;

                    //! Trial Second Moment
                    secondMoment_trial[trialTime] += std::pow(currentOrderParameter, 2.0);

                    //! Trial Mean Cluster Size
                    meanClusterSize_trial[trialTime] += model.getMeanClusterSize();
                }//* End of one step
            }//* End of network growing (one ensemble)
        } //* End of every ensembles
    } //* End of function mBFW::generate::run

    //*-------------------------------------------Save calculated variables------------------------------------------------------
    void save(){
        using namespace linearAlgebra;
        namespace fs = std::filesystem;
        fs::path p;

        //! Average and Save Order Parameter
        orderParameter /= ensembleSize;
        const std::string orderParameterPath = rootPath + "orderParameter/";
        p = orderParameterPath;
        if (!fs::exists(p)){
            fs::create_directories(p);
        }
        writeCSV(orderParameterPath+defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum, randomEngineSeed), orderParameter);

        //! Average and Save Trial Order Parameter
        orderParameter_trial /= ensembleSize;
        const std::string orderParameter_trialPath = rootPath + "orderParameter_trial/";
        p = orderParameter_trialPath;
        if (!fs::exists(p)){
            fs::create_directories(p);
        }
        writeCSV(orderParameter_trialPath+defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum, randomEngineSeed), orderParameter);

        //! Average and Save Variance
        secondMoment /= ensembleSize;
        const std::vector<double> variance = secondMoment - elementPow(orderParameter, 2.0);
        const std::string variancePath = rootPath + "variance/";
        p = variancePath;
        if (!fs::exists(p)){
            fs::create_directories(p);
        }
        writeCSV(variancePath+defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum, randomEngineSeed), variance);

        //! Average and Save Trial Variance
        secondMoment_trial /= ensembleSize;
        const std::vector<double> variance_trial = secondMoment_trial - elementPow(orderParameter_trial, 2.0);
        const std::string variance_trialPath = rootPath + "variance_trial/";
        p = variance_trialPath;
        if (!fs::exists(p)){
            fs::create_directories(p);
        }
        writeCSV(variance_trialPath+defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum, randomEngineSeed), variance_trial);

        //! Average and Save Mean Cluster Size
        meanClusterSize /= ensembleSize;
        const std::string meanClusterSizePath = rootPath + "meanClusterSize/";
        p = meanClusterSizePath;
        if (!fs::exists(p)){
            fs::create_directories(p);
        }
        writeCSV(meanClusterSizePath+defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum, randomEngineSeed), meanClusterSize);

        //! Average and Save Trial Mean Cluster Size
        meanClusterSize_trial /= ensembleSize;
        const std::string meanClusterSize_trialPath = rootPath + "meanClusterSize_trial/";
        p = meanClusterSize_trialPath;
        if (!fs::exists(p)){
            fs::create_directories(p);
        }
        writeCSV(meanClusterSize_trialPath+defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum, randomEngineSeed), meanClusterSize_trial);
    }//* End of function mBFW::generate::save
}//* End of namespace mBFW::generate