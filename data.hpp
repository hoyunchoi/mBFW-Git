#pragma once
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <set>
#include <tuple>
#include <utility>
#include <cstdio>
#include <filesystem>

#include "../library/CSV.hpp"
#include "../library/linearAlgebra.hpp"
#include "../library/dataProcess.hpp"

#include "fileName.hpp"
#include "parameters.hpp"

//? Need to seperate log binning

namespace mBFW::data{
    using namespace linearAlgebra;
    namespace fs = std::filesystem;
    const std::string rootDirectory = "../data/mBFW/";

    //* Declaration of varaibles used at mBFW::data
    const std::vector<std::string> states = {"before", "during", "after"};
    int networkSize;
    double acceptanceThreshold;
    unsigned ensembleSize;
    int coreNum;
    std::string baseName;
    int maxTime, maxTrialTime;
    bool deletion;

    //* Set parameters for average
    void setParameters(const int& t_networkSize, const double& t_acceptanceThreshold);

    //* Define additional baseDirectory
    const std::string defineAdditionalDirectory(const std::string& t_baseDirectory, const std::string& t_additional);

    //* Delete already read files
    void conditionallyDeleteFile(const std::string& t_deletionFile);

    //* Find target file at input baseDirectory
    const std::set<std::string> findTargetFileNameList(const std::string& t_directory, const std::string& t_target);

    //* Extract Ensemble Size from file name
    const unsigned extractEnsemble(const std::string& t_fileName);
    const unsigned extractTotalEnsemble(const std::set<std::string>& t_fileNameList);

    //* Extract value of repeater of standard (order parameter/time) from file name
    const double extractRepeater(const std::string& t_fileName, const std::string& t_standard, const unsigned& t_standardSize=6);

    //* Extract list of repeater of standard(order parameter/time) from file name list
    const std::set<double> extractRepeaterList(const std::set<std::string>& t_fileNameList, const std::string& t_standard);

    //* Enesemble average the input file list
    template <typename T>
    std::tuple<T, unsigned> averageFile(const std::string& t_baseDirectory, const std::set<std::string>& t_fileNameList, const T& t_dummy);

    //* Find t_a and inflection point of order parameter
    void findTa(const std::vector<double>& t_orderParameter);

    //* ----------------------------------------------------- Process for each observables --------------------------------------------------------
    //* Process (Trial)Time-X
    //! Order Parameter (Trial), Mean cluster Size (Trial), Order Parameter Variance (Trial), Inter Event Time
    void vectorAvg(const std::string& t_type);

    //* Process observable: distribution
    //! Inter Event Time Distribution (time, op), Delta Upper Bound Distribution (time, op), Age Distribution (time, op)
    template <typename T>
    void dist(const std::string& t_type, const T& t_dummy);

    //* Process observable: distribution with repeater
    //! Cluster Size Distribution (exact, time), order parameter distribution
    template <typename T>
    void repeaterDist(const std::string& t_type, const T& t_dummy);

    //* Process observable: average Y for fixed X
    void X_avgY(const std::string& t_type);

    //* Process observable: Distribution of X,Y
    void sampled_X_Y(const std::string& t_type);

    //* From sampled_X_Y, average Y for fixed X, average X for fixed Y
    void seperate(const std::string& t_type, const std::map<std::pair<int, int>, double>& t_data, const int& t_ensembleSize);

    //* Run selected observables at check list
    void run(const std::map<std::string, bool>& t_checkList);
}//* End of namespace mBFW::data

void mBFW::data::setParameters(const int& t_networkSize, const double& t_acceptanceThreshold){
    networkSize = t_networkSize;
    maxTime = t_networkSize;
    maxTrialTime = std::floor(maxTime/t_acceptanceThreshold);
    acceptanceThreshold = t_acceptanceThreshold;
    baseName = fileName::base(t_networkSize, t_acceptanceThreshold);
}

void mBFW::data::run(const std::map<std::string, bool>& t_checkList){
    if (t_checkList.at("ageDist_op")){
        for (const std::string& state : states){
            dist("ageDist_op/"+state, std::vector<double>{});
        }
    }
    if (t_checkList.at("ageDist_time")){
        for (const std::string& state : states){
            dist("ageDist_time/"+state, std::vector<double>{});
        }
    }
    if (t_checkList.at("clusterSizeDist")){
        repeaterDist("clusterSizeDist", std::map<int, double>{});
    }
    if (t_checkList.at("clusterSizeDist_exact")){
        repeaterDist("clusterSizeDist_exact", std::map<int, double>{});
    }
    if (t_checkList.at("clusterSizeDist_time")){
        repeaterDist("clusterSizeDist_time", std::map<int, double>{});
    }
    if (t_checkList.at("deltaUpperBound_deltaAcceptance")){
        for (const std::string& state : states){
            X_avgY("deltaUpperBound_deltaAcceptance/"+state);
        }
    }
    if (t_checkList.at("deltaUpperBoundDist_op")){
        for (const std::string& state : states){
            dist("deltaUpperBoundDist_op/"+state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("deltaUpperBoundDist_time")){
        for (const std::string& state : states){
            dist("deltaUpperBoundDist_time/"+state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("deltaUpperBoundDist_tot")){
        dist("deltaUpperBoundDist_tot", std::map<int, double>{});
    }
    if (t_checkList.at("interEventTime")){
        vectorAvg("interEventTime");
    }
    if (t_checkList.at("interEventTimeDist_op")){
        for (const std::string& state : states){
            dist("interEventTimeDist_op/"+state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("interEventTimeDist_time")){
        for (const std::string& state : states){
            dist("interEventTimeDist_time/"+state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("interEventTimeDist_tot")){
        dist("interEventTimeDist_tot", std::map<int, double>{});
    }
    if (t_checkList.at("interEventTime_deltaUpperBound")){
        for (const std::string& state : states){
            X_avgY("interEventTime_deltaUpperBound/"+state);
        }
    }
    if (t_checkList.at("meanClusterSize")){
        vectorAvg("meanClusterSize");
    }
    if (t_checkList.at("meanClusterSize_trial")){
        vectorAvg("meanClusterSize_trial");
    }
    if (t_checkList.at("orderParameter")){
        vectorAvg("orderParameter");
    }
    if (t_checkList.at("orderParameter_trial")){
        vectorAvg("orderParameter_trial");
    }
    if (t_checkList.at("orderParameterDist")){
        repeaterDist("orderParameterDist", std::map<double, double>{});
    }
    if (t_checkList.at("orderParameterVariance")){
        vectorAvg("orderParameterVariance");
    }
    if (t_checkList.at("orderParameterVariance_trial")){
        vectorAvg("orderParameterVariance_trial");
    }
    if (t_checkList.at("upperBound_deltaAcceptance")){
        for (const std::string& state : states){
            X_avgY("upperBound_deltaAcceptance/"+state);
        }
    }
    if (t_checkList.at("sampled_deltaUpperBound_interEventTime")){
        sampled_X_Y("sampled_deltaUpperBound_interEventTime");
    }
    if (t_checkList.at("sampled_upperBound_interEventTime")){
        sampled_X_Y("sampled_upperBound_interEventTime");
    }
    if (t_checkList.at("sampled_time_interEventTime")){
        sampled_X_Y("sampled_time_interEventTime");
    }
}

const std::string mBFW::data::defineAdditionalDirectory(const std::string& t_baseDirectory, const std::string& t_additional){
    const std::string additionalDirectory = t_baseDirectory + t_additional + "/";
    if (!fs::exists(additionalDirectory)){
        fs::create_directories(additionalDirectory);
    }
    return additionalDirectory;
}

void mBFW::data::conditionallyDeleteFile(const std::string& t_deletionFile){
    if (deletion){
        std::cout << "Deleting file " << t_deletionFile << "\n";
        CSV::deleteFile(t_deletionFile);
    }
}

const std::set<std::string> mBFW::data::findTargetFileNameList(const std::string& t_directory, const std::string& t_target){
    std::set<std::string> baseFileNameList;
    for (const auto& file : fs::directory_iterator(t_directory)){
        const std::string fileName = file.path().filename();
        if (!fileName.find(t_target)){
            baseFileNameList.emplace(fileName);
        }
    }
    return baseFileNameList;
}

const unsigned mBFW::data::extractEnsemble(const std::string& t_fileName){
    std::string temp = t_fileName.substr(t_fileName.find("E")+1);
    temp = temp.substr(0, temp.find_first_of(",-"));
    return std::stoul(temp);
}

const unsigned mBFW::data::extractTotalEnsemble(const std::set<std::string>& t_fileNameList){
    unsigned tot = 0;
    for (const std::string& fileName : t_fileNameList){
        tot += extractEnsemble(fileName);
    }
    return tot;
}

const double mBFW::data::extractRepeater(const std::string& t_fileName, const std::string& t_standard, const unsigned& t_standardSize){
    return std::stod(t_fileName.substr(t_fileName.find(t_standard) + t_standard.size(), t_standardSize));
}

const std::set<double> mBFW::data::extractRepeaterList(const std::set<std::string>& t_fileNameList, const std::string& t_standard){
    std::set<double> repeaterList;
    for (const auto& fileName : t_fileNameList){
        repeaterList.emplace(extractRepeater(fileName, t_standard));
    }
    return repeaterList;
}

template <typename T>
std::tuple<T, unsigned> mBFW::data::averageFile(const std::string& t_baseDirectory, const std::set<std::string>& t_fileNameList, const T& t_dummy){
    const unsigned totalEnsembleSize = extractTotalEnsemble(t_fileNameList);
    T average;
    for (const std::string& fileName : t_fileNameList){
        T temp;
        CSV::read(t_baseDirectory + fileName, temp);
        const double ratio = extractEnsemble(fileName)/(double)totalEnsembleSize;
        if (average.empty()){
            average = temp * ratio;
        }
        else{
            average += temp * ratio;
        }
    }
    return std::make_tuple(average, totalEnsembleSize);
}

void mBFW::data::findTa(const std::vector<double>& t_orderParameter){
    //* Smooth order parameter curve
    std::vector<double> smoothedOrderParameter(maxTime, 0.0);
    smoothedOrderParameter[0] = t_orderParameter[0];
    for (int i=1; i<maxTime-1; ++i){
        smoothedOrderParameter[i] = (t_orderParameter[i] + t_orderParameter[i-1] + t_orderParameter[i+1]) / 3.0;
    }
    smoothedOrderParameter[maxTime-1] = t_orderParameter[maxTime-1];

    //* Find maximum slope  and infelection point
    double maxSlope = 0.0;
    int inflectionIndex = 0;
    for (int i=0; i<maxTime-1; ++i){
        const double slope = smoothedOrderParameter[i+1] - smoothedOrderParameter[i];
        if (maxSlope < slope){
            maxSlope = slope;
            inflectionIndex = i;
        }
    }
    maxSlope *= networkSize;
    const std::vector<double> inflectionPoint = {(inflectionIndex+0.5)/networkSize, (smoothedOrderParameter[inflectionIndex] + smoothedOrderParameter[inflectionIndex+1])*0.5};

    //* Calculate t_a and m_a with precision
    const int L_A = (inflectionPoint[0] - inflectionPoint[1]/maxSlope) * networkSize;
    const double t_a = (double)L_A/networkSize;
    const double m_a = smoothedOrderParameter[L_A];

    //* Save the data
    std::ofstream writeFile;
    writeFile.open(rootDirectory + "points/" + fileName::base(networkSize, acceptanceThreshold) + ".txt", std::ios_base::app);
    writeFile << std::setprecision(15) << "t_a: " << t_a << "\n";
    writeFile << std::setprecision(15) << "m_a: " << m_a << "\n";
    writeFile << std::setprecision(15) << "t_inflection: " << inflectionPoint[0] << "\n";
    writeFile << std::setprecision(15) << "m_inflection: " << inflectionPoint[1] << "\n";
    writeFile.close();
}

void mBFW::data::vectorAvg(const std::string& t_type){
    //* Define Directories
    const std::string baseDirectory = defineAdditionalDirectory(rootDirectory, t_type);
    const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "average");

    //* Find target files at base directory and average directory according to system size and acceptance threshold
    const std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
    const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

    //* Check the number of files
    if (baseFileNameList.empty()){
        std::cout << "WARNING: No file at " << baseDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (additionalFileNameList.size() >= 2){
        std::cout << "WARNING: More than two files at " << additionalDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (baseFileNameList.size() == 1){
        std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
        // return;
    }

    //* Average raw files
    auto [average, totalEnsembleSize] = averageFile(baseDirectory, baseFileNameList, std::vector<double>{});
    std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
    CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

    //* Trim averaged data
    average.erase(std::remove_if(average.end()-average.size()/100, average.end(), [](const auto& value){return std::isnan(value) || value==0;}), average.end());
    for (auto& e : average){
        if (e<0){ e = 0; }
    }

    //* Write trimmed data
    std::cout << "Writing file " << additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize) << "\n";
    CSV::write(additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize), average);

    //* Find t_a for order parameter
    if (t_type == "orderParameter"){
        findTa(average);
    }

    //* Delete previous averaged and trimmed data after successfully writing
    for (const std::string& fileName : baseFileNameList){
        conditionallyDeleteFile(baseDirectory + fileName);
    }
    for (const std::string& fileName : additionalFileNameList){
        conditionallyDeleteFile(additionalDirectory + fileName);
    }
    //* void return
    return;
}

template <typename T>
void mBFW::data::dist(const std::string& t_type, const T& t_dummy){
    //* Define directories
    const std::string baseDirectory = defineAdditionalDirectory(rootDirectory, t_type);
    const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

    //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
    const std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
    const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

    //* Check the number of files
    if (baseFileNameList.empty()){
        std::cout << "WARNING: No file at " << baseDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (additionalFileNameList.size() >= 2){
        std::cout << "WARNING: More than two files at " << additionalDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (baseFileNameList.size() == 1){
        std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
        // return;
    }

    //* Average raw files and write
    auto [average, totalEnsembleSize] = averageFile(baseDirectory, baseFileNameList, t_dummy);
    average /= linearAlgebra::accumulate(average);
    std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
    CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

    //* Log Binning data and write
    const std::map<double, double> binned = dataProcess::distLogBin(average);
    std::cout << "Writing file " << additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize) << "\n";
    CSV::write(additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize), binned);

    //* Delete previous averaged and log binned data after successfully writing
    for (const std::string& fileName : baseFileNameList){
        conditionallyDeleteFile(baseDirectory + fileName);
    }
    for (const std::string fileName : additionalFileNameList){
        conditionallyDeleteFile(additionalDirectory + fileName);
    }

    //* Void return
    return;
}

template <typename T>
void mBFW::data::repeaterDist(const std::string& t_type, const T& t_dummy){
    //* Define directories
    const std::string baseDirectory = defineAdditionalDirectory(rootDirectory, t_type);
    const std::string additionalDirectory = t_type=="orderParameterDist" ? defineAdditionalDirectory(baseDirectory, "linBin") :  defineAdditionalDirectory(baseDirectory, "logBin");

    //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
    std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
    std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

    //* Decide if the cluster size distribution is accumulated by order parameter/time
    const std::string standard = (t_type.find("time") != t_type.npos || t_type == "orderParameterDist") ? "T" : "OP";

    //* Extract list of standard value = repeater from target file name list
    const std::set<double> repeaterList = extractRepeaterList(baseFileNameList, standard);

    //* Process average and log binning for every element of repeat list
    for (const double& repeater : repeaterList){
        //* Find target files at base directory and logBin directory corresponding to repeater (order parameter/time)
        std::set<std::string> baseFileNameList_repeater;
        std::set<std::string> additionalFileNameList_repeater;
        for (const std::string& fileName : baseFileNameList){
            if (fileName.find(standard + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                baseFileNameList_repeater.emplace(fileName);
            }
        }
        for (const std::string& fileName : additionalFileNameList){
            if (fileName.find(standard + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                additionalFileNameList_repeater.emplace(fileName);
            }
        }

        //* Check the number of files
        if (baseFileNameList_repeater.empty()){
            std::cout << "WARNING: No file at " << baseDirectory << ": " << baseName << "\n";
            exit(1);
        }
        else if (additionalFileNameList_repeater.size() >= 2){
            std::cout << "WARNING: More than two files at " << additionalDirectory << ": " << baseName << "\n";
            exit(1);
        }
        else if (baseFileNameList_repeater.size() == 1){
            std::cout << "Passing file " << additionalDirectory + *additionalFileNameList_repeater.begin() << "\n";
            // continue;
        }

        //* Average raw files
        auto [average, totalEnsembleSize] = averageFile(baseDirectory, baseFileNameList_repeater, t_dummy);
        average /= linearAlgebra::accumulate(average);

        //* Binning data
        const std::map<double, double> binned = t_type=="orderParameterDist" ? dataProcess::distLinBin(average) : dataProcess::distLogBin(average);

        //* Write averaged data and write log binned data
        if (standard == "OP"){
            std::cout << "Writing file " << baseDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0) << "\n";
            CSV::write(baseDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
            std::cout << "Writing file " << additionalDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater) << "\n";
            CSV::write(additionalDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);
        }
        else{
            std::cout << "Writing file " << baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0) << "\n";
            CSV::write(baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
            std::cout << "Writing file " << additionalDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater) << "\n";
            CSV::write(additionalDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);
        }

        //* Delete previous averaged and log binned data after successfully writing
        for (const std::string& fileName : baseFileNameList_repeater){
            conditionallyDeleteFile(baseDirectory + fileName);
            baseFileNameList.erase(fileName);
        }
        for (const std::string& fileName : additionalFileNameList_repeater){
            conditionallyDeleteFile(additionalDirectory + fileName);
            additionalFileNameList.erase(fileName);
        }
    }

    //* void return
    return;
}

void mBFW::data::X_avgY(const std::string& t_type){
    //* Define directories
    const std::string baseDirectory = defineAdditionalDirectory(rootDirectory, t_type);
    const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

    //* Find target files at base directory and average directory according to system size and acceptance threshold
    const std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
    const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

    //* Check the number of files
    if (baseFileNameList.empty()){
        std::cout << "WARNING: No file at " << baseDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (additionalFileNameList.size() >= 2){
        std::cout << "WARNING: More than two files at " << additionalDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (baseFileNameList.size() == 1){
        std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
        // return;
    }

    //* Average raw files and write
    auto [average, totalEnsembleSize] = averageFile(baseDirectory, baseFileNameList, std::map<int, double>{});
    std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
    CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

    //* Log Binning data
    const std::map<double, double> binned = dataProcess::avgLogBin(average);

    //* Write log binned data
    std::cout << "Writing file " << additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize) << "\n";
    CSV::write(additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize), binned);

    //* Delete previous averaged and log binned data after successfully writing
    for (const std::string& fileName : baseFileNameList){
        conditionallyDeleteFile(baseDirectory + fileName);
    }
    for (const std::string fileName : additionalFileNameList){
        conditionallyDeleteFile(additionalDirectory + fileName);
    }

    //* void return
    return;
}

void mBFW::data::sampled_X_Y(const std::string& t_type){
    //* Defind directories
    const std::string baseDirectory = defineAdditionalDirectory(rootDirectory, t_type);
    const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

    //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
    const std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
    const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

    //* Check the number of files
    if (baseFileNameList.empty()){
        std::cout << "WARNING: No file at " << baseDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (additionalFileNameList.size() >= 2){
        std::cout << "WARNING: More than two files at " << additionalDirectory << ": " << baseName << "\n";
        exit(1);
    }
    else if (baseFileNameList.size() == 1){
        std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
        // return;
    }

    //* Average raw files and write
    auto [average, totalEnsembleSize] = averageFile(baseDirectory, baseFileNameList, std::map<std::pair<int, int>, double>{});
    average /= linearAlgebra::accumulate(average);
    std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
    CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

    //* Double log binning data
    const std::map<std::pair<double, double>, double> binned = dataProcess::distLogLogBin(average);

    //* Write log binned data
    std::cout << "Writing file " << additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize) << "\n";
    CSV::write(additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize), binned);

    //* Seperate two variables
    seperate(t_type, average, totalEnsembleSize);

    //* Delete previous averaged and log binned data after successfully writing
    for (const std::string& fileName : baseFileNameList){
        conditionallyDeleteFile(baseDirectory + fileName);
    }
    for (const std::string fileName : additionalFileNameList){
        conditionallyDeleteFile(additionalDirectory + fileName);
    }

    //* void return
    return;
}


void mBFW::data::seperate(const std::string& t_type, const std::map<std::pair<int, int>, double>& t_data, const int& t_ensembleSize){
    //* Seperate two observable names
    const std::string X_Y = t_type.substr(std::string("sampled_").size());
    const std::string observable1 = X_Y.substr(0, X_Y.find("_"));
    const std::string observable2 = X_Y.substr(X_Y.find("_")+1);

    //* Define directories
    const std::string baseDirectory1 = defineAdditionalDirectory(rootDirectory, observable1 + "_" + observable2 + "_tot");
    const std::string baseDirectory2 = defineAdditionalDirectory(rootDirectory, observable2 + "_" + observable1 + "_tot");
    const std::string additionalDirectory1 = defineAdditionalDirectory(baseDirectory1, "logBin");
    const std::string additionalDirectory2 = defineAdditionalDirectory(baseDirectory2, "logBin");

    //* Find target files
    const std::set<std::string> baseFileNameList1 = findTargetFileNameList(baseDirectory1, baseName);
    const std::set<std::string> baseFileNameList2 = findTargetFileNameList(baseDirectory2, baseName);
    const std::set<std::string> additionalFileNameList1 = findTargetFileNameList(additionalDirectory1, baseName);
    const std::set<std::string> additionalFileNameList2 = findTargetFileNameList(additionalDirectory2, baseName);

    //* Check the number of files
    if (baseFileNameList1.size() >= 2){
        std::cout << "WARNING: More than two files at " << baseDirectory1 << ": " << baseName << "\n";
        exit(1);
    }
    else if (additionalFileNameList1.size() >= 2){
        std::cout << "WARNING: More than two files at " << additionalDirectory1 << ": " << baseName << "\n";
        exit(1);
    }
    else if (baseFileNameList2.size() >= 2){
        std::cout << "WARNING: More than two files at " << baseDirectory2 << ": " << baseName << "\n";
        exit(1);
    }
    else if (additionalFileNameList2.size() >= 2){
        std::cout << "WARNING: More than two files at " << additionalDirectory2 << ": " << baseName << "\n";
        exit(1);
    }

    //* Read target files if exists
    auto [average1, ensembleSize1] = averageFile(baseDirectory1, baseFileNameList1, std::map<int, double>{});
    const unsigned totalEnsembleSize1 = t_ensembleSize + ensembleSize1;
    const double oldRatio1 = ensembleSize1/(double)totalEnsembleSize1; const double newRatio1 = 1.0-oldRatio1;
    average1 *= oldRatio1;
    auto [average2, ensembleSize2] = averageFile(baseDirectory2, baseFileNameList2, std::map<int, double>{});
    const unsigned totalEnsembleSize2 = t_ensembleSize + ensembleSize2;
    const double oldRatio2 = ensembleSize2/(double)totalEnsembleSize2; const double newRatio2 = 1.0-oldRatio2;
    average2 *= oldRatio2;


    //* Get total ratio of each observables
    std::map<int, double> totalRatio1, totalRatio2;
    for (const auto& e : t_data){
        totalRatio1[e.first.first] += e.second;
        totalRatio2[e.first.second] += e.second;
    }

    //* Add value
    for (const auto& e : t_data){
        average1[e.first.first] += e.first.second * e.second/totalRatio1[e.first.first] * newRatio1;
        average2[e.first.second] += e.first.first * e.second/totalRatio2[e.first.second] * newRatio2;
    }

    //* Write averaged data
    std::cout << "Writing file " << baseDirectory1 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize1, 0) << "\n";
    CSV::write(baseDirectory1 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize1, 0), average1);
    std::cout << "Writing file " << baseDirectory2 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize2, 0) << "\n";
    CSV::write(baseDirectory2 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize2, 0), average2);

    //* Log binning data
    const std::map<double, double> binned1 = dataProcess::avgLogBin(average1);
    const std::map<double, double> binned2 = dataProcess::avgLogBin(average2);

    //* Write log binned data
    std::cout << "Writing file " << additionalDirectory1 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize1) << "\n";
    CSV::write(additionalDirectory1 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize1), binned1);
    std::cout << "Writing file " << additionalDirectory2 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize2) << "\n";
    CSV::write(additionalDirectory2 + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize2), binned2);

    //* Delete previous averaged and log binned data after successfully writing
    for (const std::string& fileName : baseFileNameList1){
        conditionallyDeleteFile(baseDirectory1 + fileName);
    }
    for (const std::string fileName : additionalFileNameList1){
        conditionallyDeleteFile(additionalDirectory1 + fileName);
    }
    for (const std::string& fileName : baseFileNameList2){
        conditionallyDeleteFile(baseDirectory2 + fileName);
    }
    for (const std::string fileName : additionalFileNameList2){
        conditionallyDeleteFile(additionalDirectory2 + fileName);
    }

    //* void return
    return;
}