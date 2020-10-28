#pragma once

#include <iostream>
#include <vector>
#include <set>

#include "../library-Git/stringFormat.hpp"

namespace mBFW{
    //! Declaration of global variables defined at mBFW namespace
    static int networkSize;
    static double acceptanceThreshold;
    static int ensembleSize;
    static int coreNum;
    static const double m_a = 0.05;
    static double m_c;
    static double t_c;
    static std::set<double> time_orderParameterDistribution;
    static std::set<double> orderParameter_clusterSizeDistribution;
    static const std::vector<std::string> states = {"before", "during"};
    const std::string rootPath = "../data/mBFW_hybrid/";

    //! Choose observables to be processed
    static bool process_orderParameter;
    static bool process_meanClusterSize;
    static bool process_secondGiant;
    static bool process_interEventTime;
    static bool process_deltaAcceptance;
    static bool process_orderParameterDistribution;
    static bool process_clusterSizeDistribution;
    static bool process_ageDistribution;
    static bool process_interEventTimeDistribution;
    static bool process_deltaUpperBoundDistribution;
    static bool process_deltaAcceptanceDistribution;
    static bool process_interEventTime_DeltaAcceptance;
    static bool process_upperBound_DeltaAcceptance;
    static bool process_deltaUpperBound_DeltaAcceptance;
    static bool process_dynamics;

    //! File name conventions
    const std::string defaultFile(const int& t_networkSize, const double& t_acceptanceThreshold){
        return "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G" + to_stringWithPrecision(t_acceptanceThreshold, 1);
    }

    const std::string defaultFileName(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const int& t_coreNum=-1, const int& t_randomEngineSeed=-1){
        const std::string fileName = "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G"+to_stringWithPrecision(t_acceptanceThreshold,1) + ",E"+std::to_string(t_ensembleSize);

        //* name for averaged or binned file
        if (t_coreNum==-1 && t_randomEngineSeed==-1){
            return fileName + ".txt";
        }

        //* name for regular file
        else if (t_randomEngineSeed == -1){
            return fileName + "-" + std::to_string(t_coreNum) + ".txt";
        }

        //* name for dynamics
        else{
            return fileName + "-" + std::to_string(t_randomEngineSeed) + ".txt";
        }
    }

    const std::string filename_time(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const double& t_time, const int& t_coreNum=-1){
        const std::string fileName = "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G"+to_stringWithPrecision(t_acceptanceThreshold,1) + ",E"+std::to_string(t_ensembleSize) + ",T"+to_stringWithPrecision(t_time,4);

        return t_coreNum==-1 ? fileName + ".txt" : fileName + "-"+std::to_string(t_coreNum)+".txt";
    }

    const std::string filename_orderParameter(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const double& t_orderParameter, const int& t_coreNum=-1){
        const std::string fileName = "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G"+to_stringWithPrecision(t_acceptanceThreshold,1) + ",E"+std::to_string(t_ensembleSize) + ",OP"+to_stringWithPrecision(t_orderParameter,4);

        return t_coreNum==-1 ? fileName + ".txt" : fileName + "-"+std::to_string(t_coreNum)+".txt";
    }

    const std::string generalFileName(const std::string& t_observable, const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const double& t_checkPoint, const int& t_coreNum=-1){
        if (t_observable == "orderParameterDistribution"){
            return filename_time(t_networkSize, t_acceptanceThreshold, t_ensembleSize, t_checkPoint, t_coreNum);
        }
        else if (t_observable == "clusterSizeDistribution"){
            return filename_orderParameter(t_networkSize, t_acceptanceThreshold, t_ensembleSize, t_checkPoint, t_coreNum);
        }
        else{
            return defaultFileName(t_networkSize, t_acceptanceThreshold, t_ensembleSize, t_coreNum);
        }
    }
} //* End of namespace mBFW