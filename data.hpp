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
#include <cstdio>
#include <filesystem>

#include "../library/CSV.hpp"
#include "../library/linearAlgebra.hpp"

#include "fileName.hpp"
#include "parameters.hpp"

namespace mBFW::data{
    using namespace linearAlgebra;
    namespace fs = std::filesystem;
    const std::string rootPath = "../data/mBFW/";

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
        std::set<std::string> baseFileNameList;
        for (const auto& file : fs::directory_iterator(t_directory)){
            const std::string fileName = file.path().filename();
            if (!fileName.find(t_target)){
                baseFileNameList.emplace(fileName);
            }
        }
        return baseFileNameList;
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
            repeaterList.emplace(extractRepeater(fileName, t_standard));
        }
        return repeaterList;
    }

    //* Integer Log Bin
    std::map<double, double> intLogBin(const std::map<int, double>& t_raw){
        //* Setup values for integer log binning
        const std::vector<double> exponentList = arange(0.0, 10.0, logBinDelta);
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
        const std::vector<double> exponentList = arange(0.0, 10.0, logBinDelta);
        const std::vector<double> min = elementPow(10.0, exponentList);
        std::vector<double> value, difference;
        value.reserve(min.size()-1); difference.reserve(min.size()-1);
        for (int i=0; i<min.size()-1; ++i){
            value.emplace_back(std::sqrt(min[i] * min[i+1]));
            difference.emplace_back(min[i+1] - min[i]);
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

    //* Double Log Bin
    std::map<double, double> doubleLogBin(const std::map<double, double>& t_raw){
        //* Setup values for double log binning
        const std::vector<double> exponentList = arange(-10.0, 0.0, logBinDelta);
        const std::vector<double> min = elementPow(10.0, exponentList);
        std::vector<double> value, difference;
        value.reserve(min.size()-1); difference.reserve(min.size()-1);
        for (int i=0; i<min.size()-1; ++i){
            value.emplace_back(std::sqrt(min[i] * min[i+1]));
            difference.emplace_back(min[i+1] - min[i]);
        }

        //* Bind the data
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

    //* Find t_a which is intercept with t-axis of tangential line of order parameter curve at the inflection point (when slope is maximum)
    void findTa(const std::vector<double>& t_orderParameter){
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
        const int l_a = (inflectionPoint[0] - inflectionPoint[1]/maxSlope) * networkSize;
        const double t_a = (double)l_a/networkSize;
        const double m_a = smoothedOrderParameter[l_a];

        //* Save the data
        std::ofstream writeFile;
        writeFile.open(rootPath + "points/" + fileName::base(networkSize, acceptanceThreshold) + ".txt", std::ios_base::app);
        writeFile << std::setprecision(15) << "t_a: " << t_a << "\n";
        writeFile << std::setprecision(15) << "m_a: " << m_a << "\n";
        writeFile << std::setprecision(15) << "t_inflection: " << inflectionPoint[0] << "\n";
        writeFile << std::setprecision(15) << "m_inflection: " << inflectionPoint[1] << "\n";
        writeFile.close();
    }

    //* Log bin inter event time with t-t_c
    void logBinIET(const std::vector<double>& t_interEventTime, const double& t_totalEnsembleSize){
        //* Change index by t-t_c
        std::map<double, double> raw;
        double t_a, m_a, t_c, m_c;
        std::tie(t_a, m_a, t_c, m_c) = mBFW::parameters::set_points(networkSize, acceptanceThreshold);
        for (int t = (int)(t_c*networkSize)+1; t<networkSize; ++t){
            raw[t/(networkSize) - t_c] = t_interEventTime[t];
        }

        //* Log Binning
        std::map<double, double> binned = doubleLogBin(raw);

        //* Check the number of files
        const std::string additionalDirectory = defineAdditionalDirectory(rootPath + "interEventTime/", "logBin");
        const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);
        if (additionalFileNameList.size() >= 2){
            std::cout << "WARNING: More than two files at " << additionalDirectory << ", " << baseName << "\n";
            exit(1);
        }

        //* Write log binned data
        std::cout << "Writing file " << additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, t_totalEnsembleSize) << "\n";
        CSV::write(additionalDirectory + fileName::NGE(networkSize, acceptanceThreshold, t_totalEnsembleSize), binned);

        //* Delete previous log binned data if writing is finished
        for (const std::string fileName : additionalFileNameList){
            conditionallyDeleteFile(additionalDirectory + fileName);
        }
    }

    //* ----------------------------------------------------- Process for each observables --------------------------------------------------------
    //* Process (Trial)Time-X
    //! Order Parameter (Trial), Mean cluster Size (Trial), Order Parameter Variance (Trial), Inter Event Time
    void time_X(const std::string& t_type){
        //* Define Directories
        const std::string baseDirectory = rootPath + t_type + "/";
        const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "average");

        //* Find target files at base directory and average directory according to system size and acceptance threshold
        const std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
        const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

        //* Check the number of files
        if (baseFileNameList.size() == 0){
            std::cout << "WARNING: No file at " << baseDirectory << ", " << baseName << "\n";
            exit(1);
        }
        else if (additionalFileNameList.size() >= 2){
            std::cout << "WARNING: More than two files at " << additionalDirectory << ", " << baseName << "\n";
            exit(1);
        }
        else if (baseFileNameList.size() == 1){
            std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
            return;
        }

        //* Extract list of ensemble size from 'targe file list' and get total ensemble size
        int totalEnsembleSize = 0;
        for (const std::string& fileName : baseFileNameList){
            totalEnsembleSize += extractEnsemble(fileName);
        }

        //* Define average vector corresponds to each observables
        std::vector<double> average;
        t_type.find("trial") != t_type.npos ? average.assign(maxTrialTime, 0.0) : average.assign(maxTime, 0.0);

        //* Read target files and average them according to weight corresponding to each ensemble size
        for (const std::string& fileName : baseFileNameList){
            std::vector<double> temp;
            CSV::read(baseDirectory + fileName, temp);
            const double ratio = extractEnsemble(fileName) / (double)totalEnsembleSize;
            average += temp * ratio;
        }

        //* Write average data
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

        //* Delete previous averaged and trimmed data after successfully writing
        for (const std::string& fileName : baseFileNameList){
            conditionallyDeleteFile(baseDirectory + fileName);
        }
        for (const std::string& fileName : additionalFileNameList){
            conditionallyDeleteFile(additionalDirectory + fileName);
        }

        //* Find t_a for order parameter
        if (t_type == "orderParameter"){
            findTa(average);
        }

        //* Do log bin for inter event time
        else if (t_type == "interEventTime"){
            logBinIET(average, totalEnsembleSize);
        }

        //* void return
        return;
    }

    //* Process Distribution distinguished by discontinuous jump
    //! Inter Event Time Distribution (time, op), Delta Upper Bound Distribution (time, op)
    void dist(const std::string& t_type){
        //* Define directories
        const std::string baseDirectory = rootPath + t_type + "/";
        const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

        //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
        const std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
        const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

        //* Check the number of files
        if (baseFileNameList.size() == 0){
            std::cout << "WARNING: No file at " << baseDirectory << ", " << baseName << "\n";
            exit(1);
        }
        else if (additionalFileNameList.size() >= 2){
            std::cout << "WARNING: More than two files at " << additionalDirectory << ", " << baseName << "\n";
            exit(1);
        }
        else if (baseFileNameList.size() == 1){
            std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
            return;
        }

        //* Find Ensemble Size of each target files and total ensemble size
        int totalEnsembleSize = 0;
        for (const std::string& fileName : baseFileNameList){
            totalEnsembleSize += extractEnsemble(fileName);
        }

        //* Read target files and average them according to weight corresponding to each ensemble size
        std::map<int, double> average;
        for (const std::string& fileName : baseFileNameList){
            std::map<int, double> temp;
            CSV::read(baseDirectory + fileName, temp);
            const double ratio = extractEnsemble(fileName) / (double)totalEnsembleSize;
            for (auto it=temp.begin(); it!= temp.end(); ++it){
                average[it->first] += it->second * ratio;
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

        //* Void return
        return;
    }

    //* Process age distribution
    //! Age Distribution (time, op)
    void ageDist(const std::string& t_type = ""){
        for (const std::string& state : states){
            //* Define directories
            const std::string baseDirectory = rootPath + "ageDist_" + t_type + "/" + state + "/";
            const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

            //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
            std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
            std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

            //* Check the number of files
            if (baseFileNameList.size() == 0){
                std::cout << "WARNING: No file at " << baseDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (additionalFileNameList.size() >= 2){
                std::cout << "WARNING: More than two files at " << additionalDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (baseFileNameList.size() == 1){
                std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
                continue;
            }

            //* Find Ensemble Size of each target files and total ensemble size
            int totalEnsembleSize = 0;
            for (const std::string& fileName : baseFileNameList){
                totalEnsembleSize += extractEnsemble(fileName);
            }

            //* Read target files and average them according to weight corresponding to each ensemble size
            std::vector<double> average(maxTime, 0.0);
            for (const std::string& fileName : baseFileNameList){
                std::vector<double> temp;
                CSV::read(baseDirectory + fileName, temp);
                const double ratio = extractEnsemble(fileName) / (double)totalEnsembleSize;
                average += temp * ratio;
            }

            //* Normalize averaged data
            average /= std::accumulate(average.begin(), average.end(), 0.0);

            //* Write averaged data
            std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
            CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

            //* Log Binning data
            std::map<double, double> binned = intLogBin(average);
            binned /= accumulate(binned);

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
        }
        //* void return
        return;
    }

    //* Process cluster size distribution
    //! Cluster Size Distribution (exact, time)
    void clustersizeDist(const std::string& t_type = ""){
        //* Define directories
        const std::string baseDirectory = rootPath + "clusterSizeDist" + t_type + "/";
        const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

        //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
        std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
        std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

        //* Decide if the cluster size distribution is accumulated by order parameter/time
        std::string standard;
        t_type.find("time") != t_type.npos ? standard = "T" : standard = "OP";

        //* Extract list of standard value = repeater from target file name list
        std::set<double> repeaterList = extractRepeaterList(baseFileNameList, standard);

        //* Process average and log binning for every element of repeat list
        for (const double& repeater : repeaterList){
            //* Find target files at base directory and logBin directory corresponding to repeater (order parameter/time)
            std::set<std::string> baseFileNameList_repeater;
            std::set<std::string> additionalFileNameList_repeater;
            int totalEnsembleSize = 0;
            for (const std::string& fileName : baseFileNameList){
                if (fileName.find(standard + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    totalEnsembleSize += extractEnsemble(fileName);
                    baseFileNameList_repeater.emplace(fileName);
                }
            }
            for (const std::string& fileName : additionalFileNameList){
                if (fileName.find(standard + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    additionalFileNameList_repeater.emplace(fileName);
                }
            }

            //* Check the number of files
            if (baseFileNameList_repeater.size() == 0){
                std::cout << "WARNING: No file at " << baseDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (additionalFileNameList_repeater.size() >= 2){
                std::cout << "WARNING: More than two files at " << additionalDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (baseFileNameList_repeater.size() == 1){
                std::cout << "Passing file " << additionalDirectory + *additionalFileNameList_repeater.begin() << "\n";
                continue;
            }

            //* Read target files and average them according to weight corresponding to each ensemble size
            std::map<int, double> average;
            for (const std::string& fileName : baseFileNameList_repeater){
                std::map<int, double> temp;
                CSV::read(baseDirectory + fileName, temp);
                const double ratio = extractEnsemble(fileName) / (double)totalEnsembleSize;
                for (auto it = temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second * ratio;
                }
            }

            //* Normalize averaged data
            average /= accumulate(average);

            //* Log Binning data
            std::map<double, double> binned = intLogBin(average);
            binned /= accumulate(binned);

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

    void opd(){
        //* Define directories
        const std::string baseDirectory = rootPath + "orderParameterDist/";
        const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "linBin");

        //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
        std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
        std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

        //* Extract list of standard value = repeater from target file name list
        std::set<double> repeaterList = extractRepeaterList(baseFileNameList, "T");

        //* Process average and log binning for every element of repeat list
        for (const double& repeater : repeaterList){
            //* Find target files at base directory and linBin directory corresponding to repeater (time)
            std::set<std::string> baseFileNameList_repeater;
            std::set<std::string> additionalFileNameList_repeater;
            int totalEnsembleSize = 0;
            for (const std::string& fileName : baseFileNameList){
                if (fileName.find("T" + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    totalEnsembleSize += extractEnsemble(fileName);
                    baseFileNameList_repeater.emplace(fileName);
                }
            }
            for (const std::string& fileName : additionalFileNameList){
                if (fileName.find("T" + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    additionalFileNameList_repeater.emplace(fileName);
                }
            }

            //* Check the number of files
            if (baseFileNameList_repeater.size() == 0){
                std::cout << "WARNING: No file at " << baseDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (additionalFileNameList_repeater.size() >= 2){
                std::cout << "WARNING: More than two files at " << additionalDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (baseFileNameList_repeater.size() == 1){
                std::cout << "Passing file " << additionalDirectory + *additionalFileNameList_repeater.begin() << "\n";
                continue;
            }

            //* Read target files and average them according to weight corresponding to each ensemble size
            std::map<double, double> average;
            for (const std::string fileName : baseFileNameList_repeater){
                std::map<double, double> temp;
                CSV::read(baseDirectory + fileName, temp);
                const double ratio = extractEnsemble(fileName) / (double)totalEnsembleSize;
                for (auto it = temp.begin(); it != temp.end(); ++it){
                    average[it->first] += it->second * ratio;
                }
            }

            //* Normalize averaged data
            average /= accumulate(average);

            //* Linear Binning data
            std::map<double, double> binned = doubleLinBin(average);
            binned /= accumulate(binned);

            //* Write averaged data and write log binned data
            std::cout << "Writing file " << baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0) << "\n";
            CSV::write(baseDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
            std::cout << "Writing file " << additionalDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater) << "\n";
            CSV::write(additionalDirectory + fileName::NGET(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);

            //* Delete previous averaged and lin binned data after successfully writing
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

    void X_deltaAcceptance(const std::string& t_type){
        for (const std::string& state : states){
            //* Define directories
            const std::string baseDirectory = rootPath + t_type + "/" + state + "/";
            const std::string additionalDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

            //* Find target files at base directory and average directory according to system size and acceptance threshold
            const std::set<std::string> baseFileNameList = findTargetFileNameList(baseDirectory, baseName);
            const std::set<std::string> additionalFileNameList = findTargetFileNameList(additionalDirectory, baseName);

            //* Check the number of files
            if (baseFileNameList.size() == 0){
                std::cout << "WARNING: No file at " << baseDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (additionalFileNameList.size() >= 2){
                std::cout << "WARNING: More than two files at " << additionalDirectory << ", " << baseName << "\n";
                exit(1);
            }
            else if (baseFileNameList.size() == 1){
                std::cout << "Passing file " << additionalDirectory + *additionalFileNameList.begin() << "\n";
                continue;
            }

            //* Extract list of ensemble size from 'targe file list' and get total ensemble size
            int totalEnsembleSize = 0;
            for (const std::string& fileName : baseFileNameList){
                totalEnsembleSize += extractEnsemble(fileName);
            }

            //* Read target files and average them according to weight corresponding to each ensemble
            std::map<int, double> average;
            for (const std::string& fileName : baseFileNameList){
                std::map<int, double> temp;
                CSV::read(baseDirectory + fileName, temp);
                const double ratio = extractEnsemble(fileName) / (double)totalEnsembleSize;
                for (auto it=temp.begin(); it!= temp.end(); ++it){
                    average[it->first] += it->second * ratio;
                }
            }

            //* Write averaged data
            std::cout << "Writing file " << baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0) << "\n";
            CSV::write(baseDirectory + fileName::NGE(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

            //* Log Binning data
            std::map<double, double> binned = intLogBin(average);

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
        }

        //* void return
        return;
    }

    void sampled_X_interEventTime(const std::string& t_type){
        //* Double log
        //* void return
        return;
    }

    //* ------------------------------------------------------------- Integration of each observables data process -----------------------------------------------------
    //* Run selected observables at check list
    void process(const std::map<std::string, bool>& t_checkList){
        if (t_checkList.at("ageDist_op")){
            ageDist("op");
        }
        if (t_checkList.at("ageDist_time")){
            ageDist("time");
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
        if (t_checkList.at("deltaUpperBound_deltaAcceptance")){
            X_deltaAcceptance("deltaUpperBound_deltaAcceptance");
        }
        if (t_checkList.at("deltaUpperBoundDist_op")){
            for (const std::string& state : states){
                dist("deltaUpperBoundDist_op/" + state);
            }
        }
        if (t_checkList.at("deltaUpperBoundDist_time")){
            for (const std::string& state : states){
                dist("deltaUpperBoundDist_time/" + state);
            }
        }
        if (t_checkList.at("deltaUpperBoundDist_tot")){
            dist("deltaUpperBoundDist_tot");
        }
        if (t_checkList.at("interEventTime")){
            time_X("interEventTime");
        }
        if (t_checkList.at("interEventTimeDist_op")){
            for (const std::string& state : states){
                dist("interEventTimeDist_op/" + state);
            }
        }
        if (t_checkList.at("interEventTimeDist_time")){
            for (const std::string& state : states){
                dist("interEventTimeDist_time/" + state);
            }
        }
        if (t_checkList.at("interEventTimeDist_tot")){
            dist("interEventTimeDist_tot");
        }
        if (t_checkList.at("interEventTime_deltaUpperBound")){
            X_deltaAcceptance("interEventTime_deltaUpperBound");
        }
        if (t_checkList.at("meanClusterSize")){
            time_X("meanClusterSize");
        }
        if (t_checkList.at("meanClusterSize_trial")){
            time_X("meanClusterSize_trial");
        }
        if (t_checkList.at("orderParameter")){
            time_X("orderParameter");
        }
        if (t_checkList.at("orderParameter_trial")){
            time_X("orderParameter_trial");
        }
        if (t_checkList.at("orderParameterDist")){
            opd();
        }
        if (t_checkList.at("orderParameterVariance")){
            time_X("orderParameterVariance");
        }
        if (t_checkList.at("orderParameterVariance_trial")){
            time_X("orderParameterVariance_trial");
        }
        if (t_checkList.at("upperBound_deltaAcceptance")){
            X_deltaAcceptance("upperBound_deltaAcceptance");
        }
        if (t_checkList.at("sampled_deltaUpperBound_interEventTime")){
            sampled_X_interEventTime("sampled_deltaUpperBound_interEventTime");
        }
        if (t_checkList.at("sampled_upperBound_interEventTime")){
            sampled_X_interEventTime("sampled_upperBound_interEventTime");
        }
        if (t_checkList.at("sampled_time_interEventTime")){
            sampled_X_interEventTime("sampled_time_interEventTime");
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