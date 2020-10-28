#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdio>
#include <filesystem>

#include "../library-Git/CSV.hpp"
#include "../library-Git/linearAlgebra.hpp"

#include "parameters.hpp"
#include "mBFW.hpp"

namespace mBFW::data{
    using namespace linearAlgebra;
    namespace fs = std::filesystem;

    //*------------------------------------------- Declaration of varaibles used at mBFW::data ------------------------------------------------------
    std::string target;
    fs::path p;
    int maxTime;
    int maxTrialTime;
    bool deletion;

    //*------------------------------------------- Set parameters for average, merge, log bin ------------------------------------------------------
    //* Set parameters for core average
    void setParameters(const int& t_networkSize, const double& t_acceptanceThreshold, const bool t_deletion){
        //! Input variables
        networkSize = t_networkSize;
        maxTime = t_networkSize;
        maxTrialTime = std::ceil(t_networkSize/t_acceptanceThreshold);
        acceptanceThreshold = t_acceptanceThreshold;
        deletion = t_deletion;
        target = defaultFile(t_networkSize, t_acceptanceThreshold);
    }

    //*------------------------------------------- functions for data process ------------------------------------------------------
    void conditionallyDeleteFile(const std::string& t_delitionFileName){
        if (deletion){
            CSV::deleteFile(t_delitionFileName);
        }
    }

    //* Extract Ensemble Size from file name
    const int extractEnsemble(const std::string t_currentFileName){
        return std::stoi(t_currentFileName.substr(target.size()+2, t_currentFileName.find_first_of(",-")));
    }

    //* Find target file at input directory
    const std::vector<std::string> findTargetFileNames(const std::string& t_directory){
        std::vector<std::string> targetFileNames;
        for (const auto& file : fs::directory_iterator(t_directory)){
            const std::string fileName = file.path().filename();
            if (!fileName.find(target)){
                targetFileNames.emplace_back(fileName);
            }
        }
        return targetFileNames;
    }

    //* Find total ensemble size
    const std::vector<int> extractEnsembleList(const std::vector<std::string>& t_targetFileNames){
        std::vector<int> ensembleList;
        for (const auto& fileName : t_targetFileNames){
            ensembleList.emplace_back(extractEnsemble(fileName));
        }
        return ensembleList;
    }

    //! Order Parameter
    void process_orderParameter(){
        //* Find target files corresponding to input System Size and Acceptance Threshold
        const std::string directory = rootPath + "orderParameter/";
        const std::vector<std::string> targetFileNames = findTargetFileNames(directory);

        //* Find Ensemble size of each files and total ensemble size
        const std::vector<int> ensembleList = extractEnsembleList(targetFileNames);
        const int totalEnsemble = std::accumulate(ensembleList.begin(), ensembleList.end(), 0);

        //* Read target files and average them according to weight
        std::vector<double> average(maxTime, 0.0);
        for (int i=0; i<targetFileNames.size(); ++i){
            std::vector<double> temp;
            CSV::read(directory + targetFileNames[i], temp);
            conditionallyDeleteFile(directory + targetFileNames[i]);
            average += temp*(double)ensembleList[i]/totalEnsemble;
        }

        //* Write averaged File
        CSV::write(directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0), average);

        //* Trim average for data plot
        const std::string average_directory = directory + "average/";
        if (!fs::exists(average_directory)){
            fs::create_directories(average_directory);
        }
        average.erase(std::remove_if(average.begin(), average.end(), [](const auto& value){return value==0;}), average.end());
        CSV::write(average_directory+defaultFileName(networkSize, acceptanceThreshold, totalEnsemble), average);
    }

    //! Order Parameter Trial
    void process_orderParameter_trial(){
        //* Find target files corresponding to input System Size and Acceptance Threshold
        const std::string directory = rootPath + "orderParameter_trial/";
        const std::vector<std::string> targetFileNames = findTargetFileNames(directory);

        //* Find Ensemble size of each files and total ensemble size
        const std::vector<int> ensembleList = extractEnsembleList(targetFileNames);
        const int totalEnsemble = std::accumulate(ensembleList.begin(), ensembleList.end(), 0);

        //* Read target files and average them according to weight
        std::vector<double> average(maxTrialTime, 0.0);
        for (int i=0; i<targetFileNames.size(); ++i){
            std::vector<double> temp;
            CSV::read(directory + targetFileNames[i], temp);
            conditionallyDeleteFile(directory + targetFileNames[i]);
            average += temp*(double)ensembleList[i]/totalEnsemble;
        }

        //* Write averaged File
        CSV::write(directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0), average);

        //* Trim average for data plot
        const std::string average_directory = directory + "average/";
        if (!fs::exists(average_directory)){
            fs::create_directories(average_directory);
        }
        average.erase(std::remove_if(average.begin(), average.end(), [](const auto& value){return value==0;}), average.end());
        CSV::write(average_directory+defaultFileName(networkSize, acceptanceThreshold, totalEnsemble), average);
    }

    //! Mean cluster Size
    void process_meanClusterSize(){
        //* Find target files corresponding to input System Size and Acceptance Threshold
        const std::string directory = rootPath + "meanClusterSize/";
        const std::vector<std::string> targetFileNames = findTargetFileNames(directory);

        //* Find Ensemble size of each files and total ensemble size
        const std::vector<int> ensembleList = extractEnsembleList(targetFileNames);
        const int totalEnsemble = std::accumulate(ensembleList.begin(), ensembleList.end(), 0);

        //* Read target files and average them according to weight
        std::vector<double> average(maxTime, 0.0);
        for (int i=0; i<targetFileNames.size(); ++i){
            std::vector<double> temp;
            CSV::read(directory + targetFileNames[i], temp);
            conditionallyDeleteFile(directory + targetFileNames[i]);
            average += temp*(double)ensembleList[i]/totalEnsemble;
        }

        //* Write averaged File
        CSV::write(directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0), average);

        //* Trim average for data plot
        const std::string average_directory = directory + "average/";
        if (!fs::exists(average_directory)){
            fs::create_directories(average_directory);
        }
        average.erase(std::remove_if(average.begin(), average.end(), [](const auto& value){return std::isnan(value) || value==0;}), average.end());
        CSV::write(average_directory+defaultFileName(networkSize, acceptanceThreshold, totalEnsemble), average);
    }

    //! Mean Cluster Size Trial
    void process_meanClusterSize_trial(){
        //* Find target files corresponding to input System Size and Acceptance Threshold
        const std::string directory = rootPath + "meanClusterSize_trial/";
        const std::vector<std::string> targetFileNames = findTargetFileNames(directory);

        //* Find Ensemble size of each files and total ensemble size
        const std::vector<int> ensembleList = extractEnsembleList(targetFileNames);
        const int totalEnsemble = std::accumulate(ensembleList.begin(), ensembleList.end(), 0);

        //* Read target files and average them according to weight
        std::vector<double> average(maxTrialTime, 0.0);
        for (int i=0; i<targetFileNames.size(); ++i){
            std::vector<double> temp;
            CSV::read(directory + targetFileNames[i], temp);
            conditionallyDeleteFile(directory + targetFileNames[i]);
            average += temp*(double)ensembleList[i]/totalEnsemble;
        }

        //* Write averaged File
        CSV::write(directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0), average);

        //* Trim average for data plot
        const std::string average_directory = directory + "average/";
        if (!fs::exists(average_directory)){
            fs::create_directories(average_directory);
        }
        average.erase(std::remove_if(average.begin(), average.end(), [](const auto& value){return std::isnan(value) || value==0;}), average.end());
        CSV::write(average_directory+defaultFileName(networkSize, acceptanceThreshold, totalEnsemble), average);
    }

    //! Order Parameter Variance
    void process_orderParameterVariance(){
        //* Find target files corresponding to input System Size and Acceptance Threshold
        const std::string directory = rootPath + "orderParameterVariance/";
        const std::vector<std::string> targetFileNames = findTargetFileNames(directory);

        //* Find Ensemble size of each files and total ensemble size
        const std::vector<int> ensembleList = extractEnsembleList(targetFileNames);
        const int totalEnsemble = std::accumulate(ensembleList.begin(), ensembleList.end(), 0);

        //* Read target files and average them according to weight
        std::vector<double> average(maxTime, 0.0);
        for (int i=0; i<targetFileNames.size(); ++i){
            std::vector<double> temp;
            CSV::read(directory + targetFileNames[i], temp);
            conditionallyDeleteFile(directory + targetFileNames[i]);
            average += temp*(double)ensembleList[i]/totalEnsemble;
        }

        //* Write averaged File
        CSV::write(directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0), average);

        //* Trim average for data plot
        const std::string average_directory = directory + "average/";
        if (!fs::exists(average_directory)){
            fs::create_directories(average_directory);
        }
        average.erase(std::remove_if(average.begin(), average.end(), [](const auto& value){value==0;}), average.end());
        for (auto& e : average){
            if (e<0){
                e = 0;
            }
        }
        CSV::write(average_directory+defaultFileName(networkSize, acceptanceThreshold, totalEnsemble), average);
    }

    //! Order Parameter Variance Trial
    void process_orderParameterVariance_trial(){
        //* Find target files corresponding to input System Size and Acceptance Threshold
        const std::string directory = rootPath + "orderParameterVariance_trial/";
        const std::vector<std::string> targetFileNames = findTargetFileNames(directory);

        //* Find Ensemble size of each files and total ensemble size
        const std::vector<int> ensembleList = extractEnsembleList(targetFileNames);
        const int totalEnsemble = std::accumulate(ensembleList.begin(), ensembleList.end(), 0);

        //* Read target files and average them according to weight
        std::vector<double> average(maxTrialTime, 0.0);
        for (int i=0; i<targetFileNames.size(); ++i){
            std::vector<double> temp;
            CSV::read(directory + targetFileNames[i], temp);
            conditionallyDeleteFile(directory + targetFileNames[i]);
            average += temp*(double)ensembleList[i]/totalEnsemble;
        }

        //* Write averaged File
        CSV::write(directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0), average);

        //* Trim average for data plot
        const std::string average_directory = directory + "average/";
        if (!fs::exists(average_directory)){
            fs::create_directories(average_directory);
        }
        average.erase(std::remove_if(average.begin(), average.end(), [](const auto& value){value==0;}), average.end());
        for (auto& e : average){
            if (e<0){
                e = 0;
            }
        }
        CSV::write(average_directory+defaultFileName(networkSize, acceptanceThreshold, totalEnsemble), average);
    }

    //* Run
    void process(){
        process_orderParameter();
        process_orderParameter_trial();
        process_meanClusterSize();
        process_meanClusterSize_trial();
        process_orderParameterVariance();
        process_orderParameterVariance_trial();
    }

    // //! average process
    // //* t_c : only for type check
    // template<typename T>
    // const std::map<T, double> average(const std::string t_directory, const T& t_check){
    //     std::map<T, double> average, temp;
    //     std::map<T, int> sampledAverage;
    //     for (int core=0; core<fileNum; ++core){
    //         const std::string readFile = t_directory + defaultFileName(networkSize, acceptanceThreshold, ensembleList[core], core);
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
    //         const std::string readFile = t_directory + generalFileName(t_observable, networkSize, acceptanceThreshold, ensembleList[core], t_checkpoint, core);
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
    //     const std::string directory = rootPath + t_observable + "/";
    //     std::vector<double> average(networkSize);
    //     std::vector<double> temp(networkSize);
    //     for (int core=0; core<fileNum; ++core){
    //         const std::string readFile = directory + defaultFileName(networkSize, acceptanceThreshold, ensembleList[core], core);
    //         CSV::read(readFile, temp);
    //         average += temp;
    //         removeFile(readFile);
    //     }
    //     average /= fileNum;

    //     const std::string writeFile1 = directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
    //     const std::string writeFile2 = directory + "average/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
    //     CSV::write(writeFile1, average);
    //     CSV::write(writeFile2, average);

    //     const std::string removeFileName = directory + "average/" + defaultFileName(networkSize, acceptanceThreshold, ensembleList[0]);
    //     removeFile(removeFileName);

    //     //* log bin w.r.t t-t_c
    //     if (t_observable == "orderParameter"){
    //         average -= m_c;
    //         std::map<double, double> raw;
    //         for (int t=int(t_c*networkSize)+1; t<networkSize; ++t){
    //             raw[(double)t/networkSize-t_c] = average[t];
    //         }
    //         const std::map<double, double> binned = logBin(raw);
    //         const std::string writeFile3 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
    //         CSV::write(writeFile3, binned);
    //         const std::string removeFileName = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleList[0]);
    //         removeFile(removeFileName);
    //     }
    // }

    // //! time vs X with log binning w.r.t t_c-t
    // void logBin_time_X(const std::string& t_observable){
    //     const std::string directory = rootPath + t_observable + "/";

    //     std::map<double, double> merge;
    //     std::map<double, int> sampledMerge;
    //     for (auto state : states){
    //         //* average
    //         const std::map<double, double> avg = average(directory + state + "/", 0.0);
    //         const std::string writeFile1 = directory + state + "/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
    //         CSV::write(writeFile1, avg);

    //         //* merge
    //         merge += avg;
    //         sampleNum(sampledMerge, avg);
    //     }
    //     merge /= sampledMerge;
    //     const std::string writeFile3 = directory + "merge/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
    //     CSV::write(writeFile3, merge);
    //     const std::string removeFileName1 = directory + "merge/" + defaultFileName(networkSize, acceptanceThreshold, ensembleList[0]);
    //     removeFile(removeFileName1);

    //     //* log bin w.r.t t_c-t
    //     merge = minus_first(t_c, merge);
    //     const std::map<double, double> binned = logBin(merge);
    //     const std::string writeFile2 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
    //     CSV::write(writeFile2, binned);
    //     const std::string removeFileName2 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleList[0]);
    //     removeFile(removeFileName2);
    // }


    // //! Check point distribution
    // //* t_c : only for type check
    // template <typename T>
    // void checkPointDistribution(const std::string& t_observable, const T& t_check){
    //     const std::string directory = rootPath + t_observable + "/";
    //     std::set<double> checkPointList;
    //     if (t_observable == "orderParameterDistribution"){
    //         checkPointList = time_orderParameterDistribution;
    //     }
    //     else if (t_observable == "clusterSizeDistribution"){
    //         checkPointList = orderParameter_clusterSizeDistribution;
    //     }

    //     for (const double& checkPoint : checkPointList){
    //         //* average
    //         const std::map<T, double> avg = averageDistribution(t_observable, directory, t_check, checkPoint);
    //         const std::string writeFile1 = directory + generalFileName(t_observable, networkSize, acceptanceThreshold, totalEnsemble, checkPoint, 0);
    //         CSV::write(writeFile1, avg);

    //         //* Linear Binning
    //         if (t_observable == "orderParameterDistribution"){
    //             const std::map<double, double> binned = linearBin(avg);
    //             const std::string writeFile2 = directory + "linearBin/" + filename_time(networkSize, acceptanceThreshold, totalEnsemble, checkPoint);
    //             CSV::write(writeFile2, avg);
    //             const std::string removeFileName2 = directory + "linearBin/" + filename_time(networkSize, acceptanceThreshold, ensembleList[0], checkPoint);
    //             removeFile(removeFileName2);
    //         }
    //         //* Log Binning
    //         else if (t_observable == "clusterSizeDistribution"){
    //             std::map<double, double> binned = logBin(avg);
    //             const double tot = accumulate(binned);
    //             binned /= tot;
    //             const std::string writeFile2 = directory + "logBin/" + filename_orderParameter(networkSize, acceptanceThreshold, totalEnsemble, checkPoint);
    //             CSV::write(writeFile2, binned);
    //             const std::string removeFileName2 = directory + "logBin/" + filename_orderParameter(networkSize, acceptanceThreshold, ensembleList[0], checkPoint);
    //             removeFile(removeFileName2);
    //         }
    //     }
    // }

    // //! Distribution of 'keys' distinguished by before and during jump
    // //* t_c : only for type check
    // template <typename T>
    // void distribution(const std::string& t_observable, const T& t_check){
    //     for (auto state : states){
    //         const std::string directory = rootPath + t_observable + "/" + state + "/";

    //         //* average
    //         const std::map<T, double> avg = averageDistribution(t_observable, directory, t_check, 0.0);
    //         const std::string writeFile1 = directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
    //         CSV::write(writeFile1, avg);

    //         //* Log binning
    //         std::map<double, double> binned = logBin(avg);
    //         const double tot = accumulate(binned);
    //         binned /= tot;
    //         const std::string writeFile2 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
    //         CSV::write(writeFile2, binned);
    //         const std::string removeFileName = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleList[0]);
    //         removeFile(removeFileName);
    //     }
    // }

    // //! X vs Delta Acceptance
    // //* t_c : only for type check
    // template <typename T>
    // void X_deltaAcceptance(const std::string&  t_observable, const T& t_check){
    //     const std::string directory = rootPath + t_observable + "/";

    //     //* average
    //     const std::map<T, double> avg = average(directory, t_check);
    //     const std::string writeFile1 = directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
    //     CSV::write(writeFile1, avg);

    //     //* Log Binning
    //     const std::map<double, double> binned = logBin(avg);
    //     const std::string writeFile2 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
    //     CSV::write(writeFile2, binned);
    //     const std::string removeFileName = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleList[0]);
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