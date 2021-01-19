#pragma once
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <set>
#include <cstdio>
#include <filesystem>

#include "../library-Git/CSV.hpp"
#include "../library-Git/linearAlgebra.hpp"

#include "fileName.hpp"

namespace mBFW::data{
    using namespace linearAlgebra;
    namespace fs = std::filesystem;
    const std::string rootPath = "../data/mBFW_hybrid/";

    //*------------------------------------------- Declaration of varaibles used at mBFW::data ------------------------------------------------------
    const std::vector<std::string> states = {"before", "during", "after"};
    int networkSize;
    double acceptanceThreshold;
    int ensembleSize;
    int coreNum;
    std::string baseName;
    int maxTime;
    int maxTrialTime;
    bool deletion;
    double logBinDelta;

    //*------------------------------------------- Set parameters for average, merge, log bin ------------------------------------------------------
    //* Set parameters for core average
    void setParameters(const int& t_networkSize, const double& t_acceptanceThreshold, const double& t_logBinDelta, const bool t_deletion){
        //! Input variables
        networkSize = t_networkSize;
        maxTime = t_networkSize;
        maxTrialTime = std::floor(maxTime/t_acceptanceThreshold);
        acceptanceThreshold = t_acceptanceThreshold;
        deletion = t_deletion;
        logBinDelta = t_logBinDelta;

        baseName = fileName::base(t_networkSize, t_acceptanceThreshold);
    }

    //*------------------------------------------- functions for data process ------------------------------------------------------

    //* Define additional baseDirectory
    const std::string defineAdditionalDirectory(const std::string& t_baseDirectory, const std::string& t_additional){
        const std::string additionalDirectory = t_baseDirectory + t_additional + "/";
        if (!fs::exists(additionalDirectory)){
            fs::create_directories(additionalDirectory);
        }
        return additionalDirectory;
    }

    //* Delete already read files
    void conditionallyDeleteFile(const std::string& t_delitionFileName){
        if (deletion){
            std::cout << "Deleting file " << t_delitionFileName << "\n";
            CSV::deleteFile(t_delitionFileName);
        }
    }

    //* Find target file at input baseDirectory
    const std::set<std::string> findTargetFileNameList(const std::string& t_directory, const std::string& t_target){
        std::set<std::string> targetFileNameList;
        for (const auto& file : fs::directory_iterator(t_directory)){
            const std::string fileName = file.path().filename();
            if (!fileName.find(t_target)){
                targetFileNameList.insert(fileName);
            }
        }
        return targetFileNameList;
    }

    //* Extract Ensemble Size from file name
    const int extractEnsemble(const std::string& t_fileName){
        std::string temp = t_fileName.substr(t_fileName.find("E")+1);
        temp = temp.substr(0, temp.find_first_of(",-"));
        return std::stoi(temp);
    }

    //* Extract value of repeater of standard (order parameter/time) from file name
    const double extractRepeater(const std::string& t_fileName, const std::string& t_standard){
        return std::stod(t_fileName.substr(t_fileName.find(t_standard) + t_standard.size(), 6));
    }

    //* Extract list of repeater of standard(order parameter/time) from file name list
    const std::set<double> extractRepeaterList(const std::set<std::string>& t_fileNameList, const std::string& t_standard){
        std::set<double> repeaterList;
        for (const auto& fileName : t_fileNameList){
            repeaterList.insert(extractRepeater(fileName, t_standard));
        }
        return repeaterList;
    }

    //* Integer Log Bin
    std::map<double, double> intLogBin(const std::map<int, double>& t_raw){
        //* Setup values for integer log binning
        const std::vector<double> exponentList = arange(0, 10, logBinDelta);
        const std::vector<double> min = elementPow(10.0, exponentList);
        std::vector<double> value, difference;
        value.reserve(min.size()-1); difference.reserve(min.size()-1);
        for (int i=0; i<min.size()-1; ++i){
            value.emplace_back(std::sqrt(min[i]*min[i+1]));
            difference.emplace_back(min[i+1]-min[i]);
        }

        //* Bin the data
        std::map<double, double> binned;
        for (auto it=t_raw.begin(); it!=t_raw.end(); ++it){
            for (int i=0; i<value.size(); ++i){
                if (it->first < min[i+1]){
                    binned[value[i]] += it->second / difference[i];
                    break;
                }
            }
        }
        return binned;
    }

    std::map<double, double> intLogBin(const std::vector<double>& t_raw){
        //* Setup values for integer log binning
        const std::vector<double> exponentList = arange(0, 10, logBinDelta);
        const std::vector<double> min = elementPow(10.0, exponentList);
        std::vector<double> value, difference;
        value.reserve(min.size()-1); difference.reserve(min.size()-1);
        for (int i=0; i<min.size()-1; ++i){
            value.emplace_back(std::sqrt(min[i]*min[i+1]));
            difference.emplace_back(min[i+1]-min[i]);
        }

        //* Bin the data
        std::map<double, double> binned;
        for (unsigned i=0; i<t_raw.size(); ++i){
            for (unsigned j=0; j<value.size(); ++j){
                if (i < min[j+1]){
                    binned[value[j]] += t_raw[i] / difference[j];
                    break;
                }
            }
        }
        return binned;
    }

    //* Double Linear Bin
    const std::map<double, double> doubleLinBin(const std::map<double, double>& t_raw){
        //* Setup values for double linear binning
        const std::vector<double> min = arange(0.0, 1.0, 5e-4);
        std::vector<double> value;
        const double difference = networkSize * 5e-4;
        for (int i=0; i<min.size()-1; ++i){
            value.emplace_back((min[i] + min[i+1])/2.0);
        }

        //* Bin the data
        std::map<double, double> binned;
        for (auto it=t_raw.begin(); it!=t_raw.end(); ++it){
            for (int i=0; i<value.size(); ++i){
                if (it->first < min[i+1]){
                    binned[value[i]] += it->second / difference;
                    break;
                }
            }
        }
        return binned;
    }

    //* ----------------------------------------------------- Process for each observables --------------------------------------------------------
    //* Process (Trial)Time-X
    //! Order Parameter (Trial), Mean cluster Size (Trial), Order Parameter Variance (Trial), Inter Event Time
    std::vector<double> time_X(const std::string& t_type){
        //* Define Directories
        const std::string baseDirectory = rootPath + t_type + "/";
        const std::string averageDirectory = defineAdditionalDirectory(baseDirectory, "average");

        //* Find target files at base directory and average directory according to system size and acceptance threshold
        const std::set<std::string> targetFileNameList = findTargetFileNameList(baseDirectory, baseName);
        const std::set<std::string> deleteFileNameList = findTargetFileNameList(averageDirectory, baseName);

        //* Break the function if only one ensemble file exists
        if (targetFileNameList.size() <= 1 && deleteFileNameList.size() <= 1){
            std::vector<double> average;
            CSV::read(averageDirectory + *deleteFileNameList.begin(), average);
            std::cout << "Passing file " << averageDirectory + *deleteFileNameList.begin() << "\n";
            return average;
        }

        //* Extract list of ensemble size from 'targe file list' and get total ensemble size
        int totalEnsembleSize = 0;
        for (const std::string& targetFileName : targetFileNameList){
            totalEnsembleSize += extractEnsemble(targetFileName);
        }

        //* Define average vector corresponds to each observables
        std::vector<double> average;
        t_type.find("trial") != t_type.npos ? average.assign(maxTrialTime, 0.0) : average.assign(maxTime, 0.0);

        //* Read target files and average them according to weight corresponding to each ensemble size
        for (const std::string& targetFileName : targetFileNameList){
            std::vector<double> temp;
            CSV::read(baseDirectory + targetFileName, temp);
            conditionallyDeleteFile(baseDirectory + targetFileName);
            average += temp * extractEnsemble(targetFileName) / (double)totalEnsembleSize;
        }

        //* Write average data
        std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
        CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

        //* Trim averaged data
        average.erase(std::remove_if(average.end()-average.size()/100, average.end(), [](const auto& value){return std::isnan(value) || value==0;}), average.end());
        for (auto& e : average){
            if (e<0){
                e = 0;
            }
        }

        //* Delete previous trimmed data
        for (const std::string& deleteFileName : deleteFileNameList){
            conditionallyDeleteFile(averageDirectory + deleteFileName);
        }

        //* Write trimmed data
        std::cout << "Writing file " << averageDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize) << "\n";
        CSV::write(averageDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize), average);

        //* Return average value of observables
        return average;
    }

    //* Process Distribution distinguished by discontinuous jump
    //! Inter Event Time Distribution (time, op), Delta Upper Bound Distribution (time, op)
    void dist(const std::string& t_type = ""){
        //* Define directories
        const std::string baseDirectory = rootPath + t_type + "/";
        const std::string logBinDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

        //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
        const std::set<std::string> targetFileNameList = findTargetFileNameList(baseDirectory, baseName);
        const std::set<std::string> deleteFileNameList = findTargetFileNameList(logBinDirectory, baseName);

        //* Break the function if only one ensemble file exists
        if (targetFileNameList.size() <= 1 && deleteFileNameList.size() <= 1){
            std::cout << "Passing file " << logBinDirectory + *deleteFileNameList.begin() << "\n";
            return;
        }

        //* Find Ensemble Size of each target files and total ensemble size
        int totalEnsembleSize = 0;
        for (const std::string& targetFileName : targetFileNameList){
            totalEnsembleSize += extractEnsemble(targetFileName);
        }

        //* Read target files and average them according to weight corresponding to each ensemble size
        std::map<int, double> average;
        for (const std::string& targetFileName : targetFileNameList){
            std::map<int, double> temp;
            CSV::read(baseDirectory + targetFileName, temp);
            conditionallyDeleteFile(baseDirectory + targetFileName);
            for (auto it=temp.begin(); it!= temp.end(); ++it){
                average[it->first] += it->second * extractEnsemble(targetFileName) / (double)totalEnsembleSize;
            }
        }

        //* Normalize averaged data
        average /= accumulate(average);

        //* Write averaged data
        std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
        CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

        //* Log Binning data
        std::map<double, double> binned = intLogBin(average);
        binned /= accumulate(binned);

        //* Delete previous log binned data
        for (const std::string fileName : deleteFileNameList){
            conditionallyDeleteFile(logBinDirectory + fileName);
        }

        //* Write log binned data
        std::cout << "Writing file " << logBinDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize) << "\n";
        CSV::write(logBinDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize), binned);

        //* void return
        return;
    }

    //* Process age distribution
    //! Age Distribution (time, op)
    void ageDist(const std::string& t_type = ""){
        for (const std::string& state : states){
            //* Define directories
            const std::string baseDirectory = rootPath + "ageDist_" + t_type + "/" + state + "/";
            const std::string logBinDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

            //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
            std::set<std::string> targetFileNameList = findTargetFileNameList(baseDirectory, baseName);
            std::set<std::string> deleteFileNameList = findTargetFileNameList(logBinDirectory, baseName);

            //* Break the function if only one ensemble file exists
            if (targetFileNameList.size() <=1 && deleteFileNameList.size() <= 1){
                std::cout << "Passing file " << logBinDirectory + *deleteFileNameList.begin() << "\n";
                continue;
            }

            //* Find Ensemble Size of each target files and total ensemble size
            int totalEnsembleSize = 0;
            for (const std::string& targetFileName : targetFileNameList){
                totalEnsembleSize += extractEnsemble(targetFileName);
            }

            //* Read target files and average them according to weight corresponding to each ensemble size
            std::vector<double> average;
            for (const std::string& targetFileName : targetFileNameList){
                std::vector<double> temp;
                CSV::read(baseDirectory + targetFileName, temp);
                conditionallyDeleteFile(baseDirectory + targetFileName);
                average += temp * extractEnsemble(targetFileName) / (double)totalEnsembleSize;
            }

            //* Normalize averaged data
            average /= std::accumulate(average.begin(), average.end(), 0.0);

            //* Write averaged data
            std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
            CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

            //* Log Binning data
            std::map<double, double> binned = intLogBin(average);
            binned /= accumulate(binned);

            //* Delete previous log binned data
            for (const std::string fileName : deleteFileNameList){
                conditionallyDeleteFile(logBinDirectory + fileName);
            }

            //* Write log binned data
            std::cout << "Writing file " << logBinDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize) << "\n";
            CSV::write(logBinDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize), binned);

            //* void return
            return;
        }
    }

    //* Process cluster size distribution
    //! Cluster Size Distribution (exact, time)
    void clustersizeDist(const std::string& t_type = ""){
        //* Define directories
        const std::string baseDirectory = rootPath + "clusterSizeDist" + t_type + "/";
        const std::string logBinDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

        //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
        std::set<std::string> targetFileNameList = findTargetFileNameList(baseDirectory, baseName);
        std::set<std::string> deleteFileNameList = findTargetFileNameList(logBinDirectory, baseName);

        //* Decide if the cluster size distribution is accumulated by order parameter/time
        std::string standard;
        t_type.find("time") != t_type.npos ? standard = "T" : standard = "OP";

        //* Extract list of standard value = repeater from target file name list
        std::set<double> repeaterList = extractRepeaterList(targetFileNameList, standard);

        //* Process average and log binning for every element of repeat list
        for (const double& repeater : repeaterList){
            //* Find target files at base directory and logBin directory corresponding to repeater (order parameter/time)
            std::set<std::string> targetFileNameList_repeater;
            std::set<std::string> deleteFileNameList_repeater;
            int totalEnsembleSize = 0;
            for (const std::string& targetFileName : targetFileNameList){
                if (targetFileName.find(standard + to_stringWithPrecision(repeater, 4)) != targetFileName.npos){
                    totalEnsembleSize += extractEnsemble(targetFileName);
                    targetFileNameList_repeater.insert(targetFileName);
                }
            }
            for (const std::string& fileName : deleteFileNameList){
                if (fileName.find(standard + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    deleteFileNameList_repeater.insert(fileName);
                }
            }
            //* Break the reapeater if only one ensemble file exists
            if (targetFileNameList_repeater.size() <= 1 && deleteFileNameList_repeater.size() <= 1){
                std::cout << "Passing file " << logBinDirectory + *deleteFileNameList_repeater.begin() << "\n";
                continue;
            }

            //* Read target files and average them according to weight corresponding to each ensemble size
            std::map<int, double> average;
            for (const std::string targetFileName : targetFileNameList_repeater){
                targetFileNameList.erase(targetFileName);
                std::map<int, double> temp;
                CSV::read(baseDirectory + targetFileName, temp);
                conditionallyDeleteFile(baseDirectory + targetFileName);
                for (auto it = temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second * extractEnsemble(targetFileName) / (double)totalEnsembleSize;
                }
            }

            //* Normalize averaged data
            average /= accumulate(average);

            //* Log Binning data
            std::map<double, double> binned = intLogBin(average);
            binned /= accumulate(binned);

            //* Delete previous log binned data
            for (const std::string& fileName : deleteFileNameList_repeater){
                conditionallyDeleteFile(logBinDirectory + fileName);
                deleteFileNameList.erase(fileName);
            }

            //* Write averaged data and write log binned data
            if (standard == "OP"){
                std::cout << "Writing file " << baseDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0) << "\n";
                CSV::write(baseDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
                std::cout << "Writing file " << logBinDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater) << "\n";
                CSV::write(logBinDirectory + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);
            }
            else{
                std::cout << "Writing file " << baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0) << "\n";
                CSV::write(baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
                std::cout << "Writing file " << logBinDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater) << "\n";
                CSV::write(logBinDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);
            }
        }

        //* void return
        return;
    }

    void opd(){
        //* Define directories
        const std::string baseDirectory = rootPath + "orderParameterDist/";
        const std::string linBinDirectory = defineAdditionalDirectory(baseDirectory, "linBin");

        //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
        std::set<std::string> targetFileNameList = findTargetFileNameList(baseDirectory, baseName);
        std::set<std::string> deleteFileNameList = findTargetFileNameList(linBinDirectory, baseName);

        //* Extract list of standard value = repeater from target file name list
        std::set<double> repeaterList = extractRepeaterList(targetFileNameList, "T");

        //* Process average and log binning for every element of repeat list
        for (const double& repeater : repeaterList){
            //* Find target files at base directory and linBin directory corresponding to repeater (time)
            std::set<std::string> targetFileNameList_repeater;
            std::set<std::string> deleteFileNameList_repeater;
            int totalEnsembleSize = 0;
            for (const std::string& targetFileName : targetFileNameList){
                if (targetFileName.find("T" + to_stringWithPrecision(repeater, 4)) != targetFileName.npos){
                    totalEnsembleSize += extractEnsemble(targetFileName);
                    targetFileNameList_repeater.insert(targetFileName);
                }
            }
            for (const std::string& fileName : deleteFileNameList){
                if (fileName.find("T" + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    deleteFileNameList_repeater.insert(fileName);
                }
            }

            //* Break the reapeater if only one ensemble file exists
            if (targetFileNameList_repeater.size() <= 1 && deleteFileNameList_repeater.size() <= 1){
                std::cout << "Passing file " << linBinDirectory + *deleteFileNameList_repeater.begin() << "\n";
                continue;
            }

            //* Read target files and average them according to weight corresponding to each ensemble size
            std::map<double, double> average;
            for (const std::string targetFileName : targetFileNameList_repeater){
                targetFileNameList.erase(targetFileName);
                std::map<double, double> temp;
                CSV::read(baseDirectory + targetFileName, temp);
                conditionallyDeleteFile(baseDirectory + targetFileName);
                for (auto it = temp.begin(); it != temp.end(); ++it){
                    average[it->first] += it->second * extractEnsemble(targetFileName) / (double)totalEnsembleSize;
                }
            }

            //* Normalize averaged data
            average /= accumulate(average);

            //* Linear Binning data
            std::map<double, double> binned = doubleLinBin(average);
            binned /= accumulate(binned);

            //* Delete previous lin binned data
            for (const std::string& fileName : deleteFileNameList_repeater){
                conditionallyDeleteFile(linBinDirectory + fileName);
                deleteFileNameList.erase(fileName);
            }

            //* Write averaged data and write log binned data
            std::cout << "Writing file " << baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0) << "\n";
            CSV::write(baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
            std::cout << "Writing file " << linBinDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater) << "\n";
            CSV::write(linBinDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);
        }

        //* void return
        return;
    }

    //* Find t_a which is intercept with t-axis of tangential line of order parameter curve at the inflection point (when slope is maximum)
    void findTa(const std::vector<double>& t_orderParameter){
        //* Smooth order parameter curve
        std::vector<double> smoothedOrderParameter(maxTime, 0.0);
        for (int i=1; i<maxTime-1; ++i){
            smoothedOrderParameter[i] = t_orderParameter[i] + t_orderParameter[i-1] + t_orderParameter[i+1];
        }
        smoothedOrderParameter /= 3;
        smoothedOrderParameter[0] = t_orderParameter[0];
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
        const int l_a = (inflectionPoint[0] - inflectionPoint[1]/maxSlope) * networkSize;
        const double t_a = (double)l_a/networkSize;
        const double m_a = smoothedOrderParameter[l_a];

        //* Print at t_a.txt
        std::ofstream writeFile;
        writeFile.open(rootPath + "t_a.txt", std::ios_base::app);
        writeFile << fileName::base(networkSize, acceptanceThreshold) << "\tinflection: " << std::setprecision(15) << inflectionPoint[0] << "," << inflectionPoint[1] << "\tt_a:" << t_a << ", m_a:" << m_a <<"\n";
        writeFile.close();
    }

    //* ------------------------------------------------------------- Integration of each observables data process -----------------------------------------------------
    //* Run selected observables at check list
    void process(const std::map<std::string, bool>& t_checkList){
        if (t_checkList.at("orderParameter")){
            std::vector<double> orderParameter = time_X("orderParameter");
            findTa(orderParameter);
        }
        if (t_checkList.at("orderParameter_trial")){
            std::vector<double> orderParameter_trial = time_X("orderParameter_trial");
        }
        if (t_checkList.at("meanClusterSize")){
            std::vector<double> meanClusterSize = time_X("meanClusterSize");
        }
        if (t_checkList.at("meanClusterSize_trial")){
            std::vector<double> meanClusterSize_trial = time_X("meanClusterSize_trial");
        }
        if (t_checkList.at("interEventTime")){
            std::vector<double> interEventTime = time_X("interEventTime");
        }
        if (t_checkList.at("orderParameterVariance")){
            std::vector<double> orderParameterVariance = time_X("orderParameterVariance");
        }
        if (t_checkList.at("orderParameterVariance_trial")){
            std::vector<double> orderParameterVariance_trial = time_X("orderParameterVariance_trial");
        }
        if (t_checkList.at("clusterSizeDist")){
            clustersizeDist();
        }
        if (t_checkList.at("clusterSizeDist_exact")){
            clustersizeDist("_exact");
        }
        if (t_checkList.at("clusterSizeDist_time")){
            clustersizeDist("_time");
        }
        if (t_checkList.at("interEventTimeDist_time")){
            for (auto state : states){
                dist("interEventTimeDist_time/" + state);
            }
        }
        if (t_checkList.at("ageDist_time")){
            for (auto state : states){
                dist("ageDist_time/" + state);
            }
        }
        if (t_checkList.at("deltaUpperBoundDist_time")){
            for (auto state : states){
                dist("deltaUpperBoundDist_time/" + state);
            }
        }
        if (t_checkList.at("orderParameterDist")){
            opd();
        }
    }

    // //! average process
    // //* t_c : only for type check
    // template<typename T>
    // const std::map<T, double> average(const std::string t_directory, const T& t_check){
    //     std::map<T, double> average, temp;
    //     std::map<T, int> sampledAverage;
    //     for (int core=0; core<fileNum; ++core){
    //         const std::string readFile = t_directory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[core], core);
    //         CSV::read(readFile, temp);
    //         average += temp;
    //         sampleNum(sampledAverage, temp);
    //         removeFile(readFile);
    //     }
    //     average /= sampledAverage;
    //     return average;
    // }

    // //! average process of distribution
    // //* for check point distribution, t_checkpoint
    // //* t_check : only for type
    // template <typename T>
    // const std::map<T, double> averageDistribution(const std::string& t_observable, const std::string& t_directory, const T& t_check, const double& t_checkpoint){
    //     std::map<T, double> average, temp;
    //     for (int core=0; core<fileNum; ++core){
    //         const std::string readFile = t_directory + generalFileName(t_observable, networkSize, acceptanceThreshold, ensembleSizeList[core], t_checkpoint, core);
    //         CSV::read(readFile, temp);
    //         average += temp;
    //         removeFile(readFile);
    //     }
    //     average /= fileNum;
    //     return average;
    // }

    // //! Linear Binning
    // template <typename T>
    // const std::map<double, double> linearBin(const std::map<T, double>& t_raw){
    //     std::map<double, double> binned;
    //     std::map<double, int> sampledBin;
    //     for (auto it=t_raw.begin(); it!=t_raw.end(); ++it){
    //         for (int i=0; i<linearBinNum; ++i){
    //             if (linear_min[i] > it->first && it->first){
    //                 binned[linear_value[i]] += it->second;
    //                 ++sampledBin[linear_value[i]];
    //             }
    //         }
    //     }
    //     binned /= sampledBin;
    //     return binned;
    // }

    // //! Log Binning
    // template <typename T>
    // const std::map<double, double> logBin(const std::map<T,double>& t_raw){
    //     //* Test whether T is double or int
    //     T test = 2;
    //     std::vector<double> min, value;
    //     if (1/test == 0.5){
    //         min = double_min;
    //         value = double_value;
    //     }
    //     else{
    //         min = int_min;
    //         value = int_value;
    //     }

    //     //* Log Binning
    //     std::map<double, double> binned;
    //     for (auto it=t_raw.begin(); it!=t_raw.end(); ++it){
    //         for (int i=0; i<logBinNum; ++i){
    //             if (min[i+1] > it->first && it->first){
    //                 binned[value[i]] += it->second/(min[i+1]-min[i]);
    //                 break;
    //             }
    //         }
    //     }
    //     return binned;
    // }

    // //! time vs X
    // void time_X(const std::string& t_observable){
    //     const std::string baseDirectory = rootPath + t_observable + "/";
    //     std::vector<double> average(networkSize);
    //     std::vector<double> temp(networkSize);
    //     for (int core=0; core<fileNum; ++core){
    //         const std::string readFile = baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[core], core);
    //         CSV::read(readFile, temp);
    //         average += temp;
    //         removeFile(readFile);
    //     }
    //     average /= fileNum;

    //     const std::string writeFile1 = baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //     const std::string writeFile2 = baseDirectory + "average/" + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile1, average);
    //     CSV::write(writeFile2, average);

    //     const std::string removeFileName = baseDirectory + "average/" + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //     removeFile(removeFileName);

    //     //* log bin w.r.t t-t_c
    //     if (t_observable == "orderParameter"){
    //         average -= m_c;
    //         std::map<double, double> raw;
    //         for (int t=int(t_c*networkSize)+1; t<networkSize; ++t){
    //             raw[(double)t/networkSize-t_c] = average[t];
    //         }
    //         const std::map<double, double> binned = logBin(raw);
    //         const std::string writeFile3 = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize);
    //         CSV::write(writeFile3, binned);
    //         const std::string removeFileName = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //         removeFile(removeFileName);
    //     }
    // }

    // //! time vs X with log binning w.r.t t_c-t
    // void logBin_time_X(const std::string& t_observable){
    //     const std::string baseDirectory = rootPath + t_observable + "/";

    //     std::map<double, double> merge;
    //     std::map<double, int> sampledMerge;
    //     for (auto state : states){
    //         //* average
    //         const std::map<double, double> avg = average(baseDirectory + state + "/", 0.0);
    //         const std::string writeFile1 = baseDirectory + state + "/" + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //         CSV::write(writeFile1, avg);

    //         //* merge
    //         merge += avg;
    //         sampleNum(sampledMerge, avg);
    //     }
    //     merge /= sampledMerge;
    //     const std::string writeFile3 = baseDirectory + "merge/" + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile3, merge);
    //     const std::string removeFileName1 = baseDirectory + "merge/" + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //     removeFile(removeFileName1);

    //     //* log bin w.r.t t_c-t
    //     merge = minus_first(t_c, merge);
    //     const std::map<double, double> binned = logBin(merge);
    //     const std::string writeFile2 = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile2, binned);
    //     const std::string removeFileName2 = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //     removeFile(removeFileName2);
    // }


    // //! Check point distribution
    // //* t_c : only for type check
    // template <typename T>
    // void checkPointDistribution(const std::string& t_observable, const T& t_check){
    //     const std::string baseDirectory = rootPath + t_observable + "/";
    //     std::set<double> checkPointList;
    //     if (t_observable == "orderParameterDistribution"){
    //         checkPointList = time_orderParameterDistribution;
    //     }
    //     else if (t_observable == "clusterSizeDistribution"){
    //         checkPointList = orderParameter_clusterSizeDistribution;
    //     }

    //     for (const double& checkPoint : checkPointList){
    //         //* average
    //         const std::map<T, double> avg = averageDistribution(t_observable, baseDirectory, t_check, checkPoint);
    //         const std::string writeFile1 = baseDirectory + generalFileName(t_observable, networkSize, acceptanceThreshold, totalEnsembleSize, checkPoint, 0);
    //         CSV::write(writeFile1, avg);

    //         //* Linear Binning
    //         if (t_observable == "orderParameterDistribution"){
    //             const std::map<double, double> binned = linearBin(avg);
    //             const std::string writeFile2 = baseDirectory + "linearBin/" + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, checkPoint);
    //             CSV::write(writeFile2, avg);
    //             const std::string removeFileName2 = baseDirectory + "linearBin/" + fileName::NGET(networkSize, acceptanceThreshold, ensembleSizeList[0], checkPoint);
    //             removeFile(removeFileName2);
    //         }
    //         //* Log Binning
    //         else if (t_observable == "clusterSizeDistribution"){
    //             std::map<double, double> binned = logBin(avg);
    //             const double tot = accumulate(binned);
    //             binned /= tot;
    //             const std::string writeFile2 = baseDirectory + "logBin/" + fileName::NGEOP(networkSize, acceptanceThreshold, totalEnsembleSize, checkPoint);
    //             CSV::write(writeFile2, binned);
    //             const std::string removeFileName2 = baseDirectory + "logBin/" + fileName::NGEOP(networkSize, acceptanceThreshold, ensembleSizeList[0], checkPoint);
    //             removeFile(removeFileName2);
    //         }
    //     }
    // }

    // //! Distribution of 'keys' distinguished by before and during jump
    // //* t_c : only for type check
    // template <typename T>
    // void distribution(const std::string& t_observable, const T& t_check){
    //     for (auto state : states){
    //         const std::string baseDirectory = rootPath + t_observable + "/" + state + "/";

    //         //* average
    //         const std::map<T, double> avg = averageDistribution(t_observable, baseDirectory, t_check, 0.0);
    //         const std::string writeFile1 = baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //         CSV::write(writeFile1, avg);

    //         //* Log binning
    //         std::map<double, double> binned = logBin(avg);
    //         const double tot = accumulate(binned);
    //         binned /= tot;
    //         const std::string writeFile2 = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize);
    //         CSV::write(writeFile2, binned);
    //         const std::string removeFileName = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //         removeFile(removeFileName);
    //     }
    // }

    // //! X vs Delta Acceptance
    // //* t_c : only for type check
    // template <typename T>
    // void X_deltaAcceptance(const std::string&  t_observable, const T& t_check){
    //     const std::string baseDirectory = rootPath + t_observable + "/";

    //     //* average
    //     const std::map<T, double> avg = average(baseDirectory, t_check);
    //     const std::string writeFile1 = baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //     CSV::write(writeFile1, avg);

    //     //* Log Binning
    //     const std::map<double, double> binned = logBin(avg);
    //     const std::string writeFile2 = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile2, binned);
    //     const std::string removeFileName = baseDirectory + "logBin/" + fileName::NGE(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //     removeFile(removeFileName);
    // }

    // //*------------------------------------------- process the data ------------------------------------------------------
    // void printProcess(const std::string& t_observable){
    //     std::cout<<"finished process : "<<t_observable<<"\n";
    // }

    // void run(){
    //     //! Order Parameter
    //     if (process_orderParameter){time_X("orderParameter"); printProcess("orderParameter");}

    //     //! Mean Cluster Size
    //     if (process_meanClusterSize){time_X("meanClusterSize"); printProcess("meanClusterSize");}

    //     //! Second giant
    //     if (process_secondGiant){time_X("secondGiant"); printProcess("secondGiant");}

    //     //! Inter Event Time
    //     if (process_interEventTime){logBin_time_X("interEventTime"); printProcess("interEventTime");}

    //     //! Delta Acceptance
    //     if (process_deltaAcceptance){logBin_time_X("deltaAcceptance"); printProcess("deltaAcceptance");}

    //     //! Order Parameter Distribution
    //     if (process_orderParameterDistribution){checkPointDistribution("orderParameterDistribution", 0.0); printProcess("orderParameterDistribution");}

    //     //! Cluster Size Distribution
    //     if (process_clusterSizeDistribution){checkPointDistribution("clusterSizeDistribution", 0); printProcess("clusterSizeDistribution");}

    //     //! Age Distribution
    //     if (process_ageDistribution){distribution("ageDistribution", 0); printProcess("orderParameter");}

    //     //! Inter Event Time Distribution
    //     if (process_interEventTimeDistribution){distribution("interEventTimeDistribution", 0); printProcess("interEventTimeDistribution");}

    //     //! Delta Upper bound Distribution
    //     if (process_deltaUpperBoundDistribution){distribution("deltaUpperBoundDistribution", 0); printProcess("deltaUpperBoundDistribution");}

    //     //! Delta Acceptance Distribution
    //     if (process_deltaAcceptanceDistribution){distribution("deltaAcceptanceDistribution", 0.0); printProcess("deltaAcceptanceDistribution");}

    //     //! Inter Event Time vs Delta Acceptance
    //     if (process_interEventTime_DeltaAcceptance){X_deltaAcceptance("interEventTime_DeltaAcceptance", 0); printProcess("interEventTime_DeltaAcceptance");}

    //     //! Upper Bound vs Delta Acceptance
    //     if (process_upperBound_DeltaAcceptance){X_deltaAcceptance("upperBound_DeltaAcceptance", 0); printProcess("upperBound_DeltaAcceptance");}

    //     //! Delta Upper Bound vs Delta Acceptance
    //     if (process_deltaUpperBound_DeltaAcceptance){X_deltaAcceptance("deltaUpperBound_DeltaAcceptance", 0); printProcess("deltaUpperBound_DeltaAcceptance");}
    // }
}//* End of namespace mBFW::data