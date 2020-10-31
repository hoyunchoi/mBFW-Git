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
    double logBinDelta;

    //*------------------------------------------- Set parameters for average, merge, log bin ------------------------------------------------------
    //* Set parameters for core average
    void setParameters(const int& t_networkSize, const double& t_acceptanceThreshold, const double& t_logBinDelta, const bool t_deletion){
        //! Input variables
        networkSize = t_networkSize;
        maxTime = t_networkSize;
        maxTrialTime = std::ceil(t_networkSize/t_acceptanceThreshold);
        acceptanceThreshold = t_acceptanceThreshold;
        deletion = t_deletion;
        target = defaultFile(t_networkSize, t_acceptanceThreshold);
        logBinDelta = t_logBinDelta;

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
    const std::vector<std::string> findTargetFileNameList(const std::string& t_directory){
        std::vector<std::string> targetFileNameList;
        for (const auto& file : fs::directory_iterator(t_directory)){
            const std::string fileName = file.path().filename();
            if (!fileName.find(target)){
                targetFileNameList.emplace_back(fileName);
            }
        }
        return targetFileNameList;
    }

    //* Extract Ensemble Size from file name
    const int extractEnsemble(std::string t_fileName){
        const int index = t_fileName.find(target);
        t_fileName = t_fileName.substr(index + target.size()+2);
        return std::stoi(t_fileName.substr(0, t_fileName.find_first_of(",-")));
    }

    //* Extract Ensemble List from target file name list
    const std::vector<int> extractEnsembleList(const std::vector<std::string>& t_fileNameList){
        std::vector<int> ensembleSizeList;
        for (const auto& fileName : t_fileNameList){
            ensembleSizeList.emplace_back(extractEnsemble(fileName));
        }
        return ensembleSizeList;
    }

    //* Extract value of repeater of standard (order parameter/time) from file name
    const double extractRepeater(const std::string& t_fileName, const std::string& t_standard){
        return std::stod(t_fileName.substr(t_fileName.find(t_standard) + t_standard.size(), 6));
    }

    //* Extract list of repeater of standard(order parameter/time) from file name list
    const std::set<double> extractRepeaterList(const std::vector<std::string>& t_fileNameList, const std::string& t_standard){
        std::set<double> repeaterList;
        for (const auto& fileName : t_fileNameList){
            repeaterList.insert(extractRepeater(fileName, t_standard));
        }
        return repeaterList;
    }

    //* Integer Log Bin
    const std::map<double, double> intLogBin(const std::map<int, double>& t_raw){
        //* Setup values for integer log binning
        const std::vector<double> exponentList = arange(0, 10, logBinDelta);
        const std::vector<double> min = elementPow(10.0, exponentList);
        std::vector<double> value, difference;
        for (int i=0; i<min.size()-1; ++i){
            value.emplace_back(std::sqrt(min[i]*min[i+1]));
            difference.emplace_back(min[i+1]-min[i]);
        }

        //* Bin the data
        std::map<double, double> binned;
        for (auto it=t_raw.begin(); it!=t_raw.end(); ++it){
            for (int i=0; i<min.size()-1; ++i){
                if (it->first < min[i+1]){
                    binned[value[i]] += it->second/difference[i];
                    break;
                }
            }
        }

        return binned;
    }

    //* Process Time-X (Trial)
    //! Order Parameter (Trial), Mean cluster Size (Trial), Order Parameter Variance (Trial)
    void time_X(const std::string& t_type){
        //* Define Directories
        const std::string baseDirectory = rootPath + t_type + "/";
        const std::string averageDirectory = defineAdditionalDirectory(baseDirectory, "average");

        //* Find target files at base directory and average directory according to system size and acceptance threshold
        const std::vector<std::string> targetFileNameList = findTargetFileNameList(baseDirectory);
        const std::vector<std::string> deleteFileNameList = findTargetFileNameList(averageDirectory);

        //* Extract list of ensemble size from 'targe file list' and get total ensemble size
        const std::vector<int> ensembleSizeList = extractEnsembleList(targetFileNameList);
        const int totalEnsembleSize = std::accumulate(ensembleSizeList.begin(), ensembleSizeList.end(), 0);

        //* Define average vector
        std::vector<double> average;
        t_type.find("trial") != t_type.npos ? average.assign(maxTrialTime, 0.0) : average.assign(maxTime, 0.0);


        //* Read target files and average them according to weight corresponding to each ensemble size
        for (int i=0; i<targetFileNameList.size(); ++i){
            std::vector<double> temp;
            CSV::read(baseDirectory + targetFileNameList[i], temp);
            conditionallyDeleteFile(baseDirectory + targetFileNameList[i]);
            average += temp*(double)ensembleSizeList[i]/totalEnsembleSize;
        }

        //* Write averaged file into base diretory
        CSV::write(baseDirectory + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize, 0), average);

        //* Trim averaged data
        average.erase(std::remove_if(average.begin(), average.end(), [](const auto& value){return std::isnan(value) || value==0;}), average.end());
        for (auto& e : average){
            if (e<0){
                e = 0;
            }
        }

        //* Delete previous trimmed data and write new trimmed data
        for (const std::string& deleteFileName : deleteFileNameList){
            conditionallyDeleteFile(averageDirectory + deleteFileName);
        }
        CSV::write(averageDirectory + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize), average);
    }

    //* Process cluster size distribution
    //! Cluster Size Distribution (exact, time)
    void clustersizeDist(const std::string& t_type = ""){
        //* Define directories
        const std::string baseDirectory = rootPath + "clusterSizeDist" + t_type + "/";
        const std::string logBinDirectory = defineAdditionalDirectory(baseDirectory, "logBin");

        //* Find target files at base directory and logBin directory corresponding to input system size and acceptance threshold
        const std::vector<std::string> targetFileNameList = findTargetFileNameList(baseDirectory);
        const std::vector<std::string> deleteFileNameList = findTargetFileNameList(logBinDirectory);

        //* Decide if the cluster size distribution is accumulated by order parameter/time
        std::string standard;
        t_type.find("time") != t_type.npos ? standard = "T" : standard = "OP";

        //* Extract list of standard value = repeater from target file name list
        std::set<double> repeaterList = extractRepeaterList(targetFileNameList, standard);

        //* Process average and log binning for every element of repeat list
        for (const double& repeater : repeaterList){
            //* Find target files at base directory and logBin directory corresponding to repeater (order parameter/time)
            std::vector<std::string> targetFileNameList_repeater;
            std::vector<std::string> deleteFileNameList_repeater;
            for (const std::string& fileName : targetFileNameList){
                if (fileName.find(standard + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    targetFileNameList_repeater.emplace_back(fileName);
                }
            }
            for (const std::string& fileName : deleteFileNameList){
                if (fileName.find(standard + to_stringWithPrecision(repeater, 4)) != fileName.npos){
                    deleteFileNameList_repeater.emplace_back(fileName);
                }
            }

            //* Find Ensemble Size of each target files and total ensemble size
            const std::vector<int> ensembleSizeList = extractEnsembleList(targetFileNameList_repeater);
            const int totalEnsembleSize = std::accumulate(ensembleSizeList.begin(), ensembleSizeList.end(), 0);

            //* Break the reapeater if only one ensemble file exists
            if (targetFileNameList_repeater.size() == 1){
                break;
            }

            //* Read target files and average them according to weight corresponding to each ensemble size
            std::map<int, double> average;
            for (int i=0; i<targetFileNameList_repeater.size(); ++i){
                std::map<int, double> temp;
                CSV::read(baseDirectory + targetFileNameList_repeater[i], temp);
                conditionallyDeleteFile(baseDirectory + targetFileNameList_repeater[i]);
                for (auto it=temp.begin(); it!= temp.end(); ++it){
                    average[it->first] += it->second*ensembleSizeList[i]/totalEnsembleSize;
                }
            }

            //* Normalize averaged data and Log Binned data
            average /= accumulate(average);
            std::map<double, double> binned = intLogBin(average);
            const double tot = accumulate(binned);
            binned /= tot;

            //* Delete previous log binned data
            for (const std::string fileName : deleteFileNameList_repeater){
                conditionallyDeleteFile(logBinDirectory + fileName);
            }

            //* Write averaged data and write log binned data
            if (standard == "OP"){
                CSV::write(baseDirectory + filename_orderParameter(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
                CSV::write(logBinDirectory + filename_orderParameter(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);
            }
            else{
                CSV::write(baseDirectory + filename_time(networkSize, acceptanceThreshold, totalEnsembleSize, repeater, 0), average);
                CSV::write(logBinDirectory + filename_time(networkSize, acceptanceThreshold, totalEnsembleSize, repeater), binned);
            }
        }
    }

    //* Run Every observables
    void process(){
        time_X("orderParameter");
        time_X("orderParameter_trial");
        time_X("meanClusterSize");
        time_X("meanClusterSize_trial");
        time_X("orderParameterVariance");
        time_X("orderParameterVariance_trial");
        clustersizeDist();
        clustersizeDist("_exact");
        clustersizeDist("_time");

    }

    // //! average process
    // //* t_c : only for type check
    // template<typename T>
    // const std::map<T, double> average(const std::string t_directory, const T& t_check){
    //     std::map<T, double> average, temp;
    //     std::map<T, int> sampledAverage;
    //     for (int core=0; core<fileNum; ++core){
    //         const std::string readFile = t_directory + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[core], core);
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
    //         const std::string readFile = baseDirectory + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[core], core);
    //         CSV::read(readFile, temp);
    //         average += temp;
    //         removeFile(readFile);
    //     }
    //     average /= fileNum;

    //     const std::string writeFile1 = baseDirectory + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //     const std::string writeFile2 = baseDirectory + "average/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile1, average);
    //     CSV::write(writeFile2, average);

    //     const std::string removeFileName = baseDirectory + "average/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //     removeFile(removeFileName);

    //     //* log bin w.r.t t-t_c
    //     if (t_observable == "orderParameter"){
    //         average -= m_c;
    //         std::map<double, double> raw;
    //         for (int t=int(t_c*networkSize)+1; t<networkSize; ++t){
    //             raw[(double)t/networkSize-t_c] = average[t];
    //         }
    //         const std::map<double, double> binned = logBin(raw);
    //         const std::string writeFile3 = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize);
    //         CSV::write(writeFile3, binned);
    //         const std::string removeFileName = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[0]);
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
    //         const std::string writeFile1 = baseDirectory + state + "/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //         CSV::write(writeFile1, avg);

    //         //* merge
    //         merge += avg;
    //         sampleNum(sampledMerge, avg);
    //     }
    //     merge /= sampledMerge;
    //     const std::string writeFile3 = baseDirectory + "merge/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile3, merge);
    //     const std::string removeFileName1 = baseDirectory + "merge/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[0]);
    //     removeFile(removeFileName1);

    //     //* log bin w.r.t t_c-t
    //     merge = minus_first(t_c, merge);
    //     const std::map<double, double> binned = logBin(merge);
    //     const std::string writeFile2 = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile2, binned);
    //     const std::string removeFileName2 = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[0]);
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
    //             const std::string writeFile2 = baseDirectory + "linearBin/" + filename_time(networkSize, acceptanceThreshold, totalEnsembleSize, checkPoint);
    //             CSV::write(writeFile2, avg);
    //             const std::string removeFileName2 = baseDirectory + "linearBin/" + filename_time(networkSize, acceptanceThreshold, ensembleSizeList[0], checkPoint);
    //             removeFile(removeFileName2);
    //         }
    //         //* Log Binning
    //         else if (t_observable == "clusterSizeDistribution"){
    //             std::map<double, double> binned = logBin(avg);
    //             const double tot = accumulate(binned);
    //             binned /= tot;
    //             const std::string writeFile2 = baseDirectory + "logBin/" + filename_orderParameter(networkSize, acceptanceThreshold, totalEnsembleSize, checkPoint);
    //             CSV::write(writeFile2, binned);
    //             const std::string removeFileName2 = baseDirectory + "logBin/" + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSizeList[0], checkPoint);
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
    //         const std::string writeFile1 = baseDirectory + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //         CSV::write(writeFile1, avg);

    //         //* Log binning
    //         std::map<double, double> binned = logBin(avg);
    //         const double tot = accumulate(binned);
    //         binned /= tot;
    //         const std::string writeFile2 = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize);
    //         CSV::write(writeFile2, binned);
    //         const std::string removeFileName = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[0]);
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
    //     const std::string writeFile1 = baseDirectory + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize, 0);
    //     CSV::write(writeFile1, avg);

    //     //* Log Binning
    //     const std::map<double, double> binned = logBin(avg);
    //     const std::string writeFile2 = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsembleSize);
    //     CSV::write(writeFile2, binned);
    //     const std::string removeFileName = baseDirectory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSizeList[0]);
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