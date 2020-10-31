# pragma once

#include <iostream>
#include <fstream>
#include <set>

#include "fileName.hpp"

namespace mBFW::parameters{
    const std::string rootPath = "../data/mBFW_hybrid/";

    //* Set parameter for orderParameter_clusterSizeDist
    std::set<double> set_orderParameter_clusterSizeDist(const int& t_networkSize, const double& t_acceptanceThreshold){
        std::set<double> orderParameter_clusterSizeDist;
        if (t_acceptanceThreshold == 0.2){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.992, 0.994, 0.996, 0.998};
        }
        else if (t_acceptanceThreshold == 0.3){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.955, 0.96, 0.965, 0.97};
        }
        else if (t_acceptanceThreshold == 0.4){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.905, 0.91, 0.915, 0.92, 0.95};
        }
        else if (t_acceptanceThreshold == 0.5){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9};
        }
        else if (t_acceptanceThreshold == 0.6){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.9};
        }
        else if (t_acceptanceThreshold == 0.7){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.67, 0.69, 0.7, 0.71, 0.73, 0.8, 0.9};
        }
        else if (t_acceptanceThreshold == 0.8){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.7, 0.8, 0.9};
        }
        else if (t_acceptanceThreshold == 0.9){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.5, 0.6, 0.7, 0.8, 0.9};
        }
        else if (t_acceptanceThreshold == 1.0){
            orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
        }
        else{
            std::cout << "Not defined acceptance rate\n";
        }
        return orderParameter_clusterSizeDist;
    }

    //* Set parameter for time_clusterSizeDist
    std::set<double> set_time_clusterSizeDist(const int& t_networkSize, const double& t_acceptanceThreshold){
        std::set<double> time_clusterSizeDist;
        if (t_acceptanceThreshold == 0.2){
            time_clusterSizeDist = {0.98, 0.985, 0.99, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999};
        }
        else if (t_acceptanceThreshold == 0.3){
            time_clusterSizeDist = {0.95, 0.96, 0.97, 0.98, 0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.989, 0.99, 0.995};
        }
        else if (t_acceptanceThreshold == 0.4){
            time_clusterSizeDist = {0.92, 0.94, 0.96, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969, 0.97, 0.98, 0.99};
        }
        else if (t_acceptanceThreshold == 0.5){
            time_clusterSizeDist = {0.9, 0.91, 0.92, 0.93, 0.931, 0.932, 0.933, 0.934, 0.935, 0.936, 0.937, 0.938, 0.939, 0.94, 0.95, 0.96};
        }
        else if (t_acceptanceThreshold == 0.6){
            time_clusterSizeDist = {0.84, 0.86, 0.88, 0.89, 0.896, 0.898, 0.90, 0.902, 0.904, 0.905, 0.906, 0.908, 0.91, 0.92, 0.94};
        }
        else if (t_acceptanceThreshold == 0.7){
            time_clusterSizeDist = {0.8, 0.82, 0.83, 0.84, 0.846, 0.848, 0.85, 0.852, 0.854, 0.856, 0.858, 0.86, 0.87,};
        }
        else if (t_acceptanceThreshold == 0.8){
            time_clusterSizeDist = {0.74, 0.76, 0.78, 0.79, 0.792, 0.794, 0.796, 0.798, 0.80, 0.802, 0.804, 0.806, 0.808, 0.81, 0.83};
        }
        else if (t_acceptanceThreshold == 0.9){
            time_clusterSizeDist = {0.66, 0.68, 0.69, 0.7, 0.71, 0.712, 0.714, 0.716, 0.718, 0.72, 0.722, 0.724, 0.726, 0.728, 0.73, 0.74};
        }
        else if (t_acceptanceThreshold == 1.0){
            time_clusterSizeDist = {0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55};
        }
        else{
            std::cout << "Not defined acceptance rate\n";
        }
        return time_clusterSizeDist;
    }

    //* Read t_a from "orderParameter/inflection.txt" and set t_a
    double set_ta(const int& t_networkSize, const double& t_acceptanceThreshold){
        double t_a = 0.0;
        std::ifstream inflectionFile(rootPath + "orderParameter/inflection.txt");
        std::string line;
        while (getline(inflectionFile, line)){
            //* Find line for input network size and acceptance threshold
            if (line.find(fileName::base(t_networkSize, t_acceptanceThreshold)) != line.npos){
                line = line.substr(line.find("t_a")+5);
                t_a = std::stod(line);
            }
        }
        if (!t_a){
            std::cout << "Can't find t_a\n";
        }
        return t_a;
    }




    // std::tuple<double, double, std::set<double>, std::set<double>> pre_defined(const int&t_networkSize, const double& t_acceptanceThreshold){
    //     if (t_acceptanceThreshold == 0.2){
    //         m_c = 0.989; t_c = 0.995;
    //         time_orderParameterDistribution = {0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985,0.999,0.9995};
    //         orderParameter_clusterSizeDistribution = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //     }
    //     else if (t_acceptanceThreshold == 0.3){
    //         m_c = 0.960; t_c = 0.983;
    //     }
    //     else if (t_acceptanceThreshold == 0.4){
    //         m_c = 0.910; t_c = 0.964;
    //     }
    //     else if (t_acceptanceThreshold == 0.5){
    //         m_c = 0.860; t_c = 0.937;
    //         time_orderParameterDistribution = {0.9320, 0.9370, 0.9482};
    //         orderParameter_clusterSizeDistribution = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.8254, 0.86, 0.9};
    //         // orderParameter_clusterSizeDistribution = {0.2};
    //     }
    //     else if (t_acceptanceThreshold == 0.6){
    //         m_c = 0.784; t_c = 0.901;
    //     }
    //     else if (t_acceptanceThreshold == 0.7){
    //         m_c = 0.693; t_c = 0.856;
    //     }
    //     else if (t_acceptanceThreshold == 0.8){
    //         m_c = 0.580; t_c = 0.799;
    //     }
    //     else if (t_acceptanceThreshold == 0.9){
    //         m_c = 0.415; t_c = 0.718;
    //     }
    //     else if (t_acceptanceThreshold == 1.0){
    //         m_c = 0.000; t_c = 0.500;
    //     }
    //     return std::make_tuple(m_c, t_c, time_orderParameterDistribution, orderParameter_clusterSizeDistribution);
    // }//* End of function mBFW::parameters::pre_defined



    // std::tuple<std::vector<double>, std::vector<double>, double, double> getParameters(const int &t_networkSize, const double &t_g){
    //     std::vector<int> N={10000,20000,40000,80000,160000,320000,640000,1280000,2560000,5120000,10240000};
    //     std::map<int,std::vector<double>> timeOfOrderParameterDistribution;
    //     std::map<int, std::vector<double>> orderParameterOfClusterSizeDistribution;
    //     double m_c, t_c;

    //     if (t_g==0.5){
    //         m_c=0.86;
    //         t_c = 0.937;

    //         timeOfOrderParameterDistribution[10000]=std::vector<double>{0.9536,0.9537,0.9538,0.9539,0.954,0.9541,0.9542,0.9543,0.9544,0.9545};
    //         timeOfOrderParameterDistribution[20000]=std::vector<double>{0.9501,0.9502,0.9503,0.9504};
    //         timeOfOrderParameterDistribution[40000]=std::vector<double>{0.9481,0.9482,0.9483,0.9484};
    //         timeOfOrderParameterDistribution[80000]=std::vector<double>{0.9466,0.9467,0.9468,0.9469};
    //         timeOfOrderParameterDistribution[160000]=std::vector<double>{0.9461,0.9462,0.9463,0.9464};
    //         timeOfOrderParameterDistribution[320000]=std::vector<double>{0.9451,0.9452,0.9453,0.9454};
    //         timeOfOrderParameterDistribution[640000]=std::vector<double>{0.9446,0.9447,0.9448,0.9449};
    //         timeOfOrderParameterDistribution[1280000]=std::vector<double>{0.9446,0.9447,0.9448,0.9449};
    //         timeOfOrderParameterDistribution[2560000]=std::vector<double>{0.9441,0.9442,0.9443,0.9444};
    //         timeOfOrderParameterDistribution[5120000]=std::vector<double>{0.9441,0.9442,0.9443,0.9444};
    //         timeOfOrderParameterDistribution[10240000]=std::vector<double>{0.9441,0.9442,0.9443,0.9444};

    //         orderParameterOfClusterSizeDistribution[10000]=std::vector<double>{0.8502,0.8504,0.8506,0.8508,0.851,0.8512,0.8514,0.8516,0.8518,0.852,0.8522,0.8524,0.8526,0.8528,0.853,0.8532,0.8534,0.8536,0.8538,0.854,0.8542,0.8544,0.8546};
    //         orderParameterOfClusterSizeDistribution[20000]=std::vector<double>{0.8502,0.8504,0.8506,0.8508,0.851,0.8512,0.8514,0.8516,0.8518,0.852,0.8522,0.8524,0.8526,0.8528,0.853,0.8532,0.8534,0.8536,0.8538,0.854,0.8542,0.8544,0.8546};
    //         orderParameterOfClusterSizeDistribution[40000]=std::vector<double>{0.8502,0.8504,0.8506,0.8508,0.851,0.8512,0.8514,0.8516,0.8518,0.852,0.8522,0.8524,0.8526,0.8528,0.853,0.8532,0.8534,0.8536,0.8538,0.854,0.8542,0.8544,0.8546,0.8548};
    //         orderParameterOfClusterSizeDistribution[80000]=std::vector<double>{0.8452,0.8454,0.8456,0.8458,0.846,0.8462,0.8464,0.8466,0.8468,0.847,0.8472,0.8474,0.8476,0.8478,0.848,0.8482,0.8484,0.8486,0.8488,0.849,0.8492,0.8494,0.8496,0.8498};
    //         orderParameterOfClusterSizeDistribution[160000]=std::vector<double>{0.842,0.8422,0.8424,0.8426,0.8428,0.843,0.8432,0.8434,0.8436,0.8438,0.844,0.8442,0.8444,0.8446,0.8448,0.8452,0.8454,0.8456,0.8458,0.846,0.8462,0.8464,0.8466,0.8468,0.847,0.8472,0.8474,0.8476,0.8478,0.848};
    //         orderParameterOfClusterSizeDistribution[320000]=std::vector<double>{0.8352,0.8354,0.8356,0.8358,0.836,0.8362,0.8364,0.8366,0.8368,0.837,0.8372,0.8374,0.8376,0.8378,0.838,0.8382,0.8384,0.8386,0.8388,0.839,0.8392,0.8394,0.8396,0.8398};
    //         orderParameterOfClusterSizeDistribution[640000]=std::vector<double>{0.827,0.8272,0.8274,0.8276,0.8278,0.828,0.8282,0.8284,0.8286,0.8288,0.829,0.8292,0.8294,0.8296,0.8298,0.8302,0.8304,0.8306,0.8308,0.831,0.8312,0.8314,0.8316,0.8318,0.832,0.8322,0.8324,0.8326,0.8328,0.833};
    //         orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.15,0.2,0.3,0.4,0.8254};
    //         orderParameterOfClusterSizeDistribution[2560000]=std::vector<double>{0.8252,0.8254,0.8256,0.8258,0.826,0.8262,0.8264,0.8266,0.8268,0.827,0.8272,0.8274,0.8276,0.8278,0.828,0.8282,0.8284,0.8286,0.8288,0.829,0.8292,0.8294,0.8296,0.8298};
    //         orderParameterOfClusterSizeDistribution[5120000]=std::vector<double>{0.8152,0.8154,0.8156,0.8158,0.816,0.8162,0.8164,0.8166,0.8168,0.817,0.8172,0.8174,0.8176,0.8178,0.818,0.8182,0.8184,0.8186,0.8188,0.819,0.8192,0.8194,0.8196,0.8198};
    //         orderParameterOfClusterSizeDistribution[10240000]=std::vector<double>{0.8152,0.8154,0.8156,0.8158,0.816,0.8162,0.8164,0.8166,0.8168,0.817,0.8172,0.8174,0.8176,0.8178,0.818,0.8182,0.8184,0.8186,0.8188,0.819,0.8192,0.8194,0.8196,0.8198};
    //     }
    //     else if (t_g==0.2){
    //         m_c=0.989;
    //         t_c = 0.995;

    //         timeOfOrderParameterDistribution[10000]=std::vector<double> {0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985,0.999,0.9995};
    //         timeOfOrderParameterDistribution[20000]=std::vector<double> {0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985,0.999,0.9995};
    //         timeOfOrderParameterDistribution[40000]=std::vector<double> {0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985,0.999,0.9995};
    //         timeOfOrderParameterDistribution[80000]=std::vector<double> {0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985,0.999};
    //         timeOfOrderParameterDistribution[160000]=std::vector<double> {0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985,0.999};
    //         timeOfOrderParameterDistribution[320000]=std::vector<double> {0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985,0.999};
    //         timeOfOrderParameterDistribution[640000]=std::vector<double> {0.992,0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985};
    //         timeOfOrderParameterDistribution[1280000]=std::vector<double> {0.992,0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985};
    //         timeOfOrderParameterDistribution[2560000]=std::vector<double> {0.992,0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998,0.9985};
    //         timeOfOrderParameterDistribution[5120000]=std::vector<double> {0.9915,0.992,0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998};
    //         timeOfOrderParameterDistribution[10240000]=std::vector<double> {0.9915,0.992,0.9925,0.993,0.9935,0.994,0.9945,0.995,0.9955,0.996,0.9965,0.997,0.9975,0.998};

    //         orderParameterOfClusterSizeDistribution[10000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[20000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[40000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[80000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[160000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[320000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[640000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[2560000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[5120000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[10240000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.95,0.99};

    //     }
    //     else if (t_g==0.3){
    //         m_c=0.96;
    //         t_c = 0.983;

    //         timeOfOrderParameterDistribution[10000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[20000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[40000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[80000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[160000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[320000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[640000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[1280000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[2560000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[5120000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         timeOfOrderParameterDistribution[10240000]=std::vector<double>{0.96,0.97,0.975,0.98,0.981,0.982,0.983,0.984,0.985,0.986,0.987,0.988,0.989,0.99,0.991,0.992,0.993,0.994,0.995};
    //         orderParameterOfClusterSizeDistribution[10000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[20000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[40000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[80000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[160000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[320000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[640000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         // orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.95,0.99};
    //         orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[2560000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[5120000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         orderParameterOfClusterSizeDistribution[10240000]=std::vector<double>{0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};

    //     }
    //     else if (t_g==0.4){
    //         m_c=0.91;
    //         t_c = 0.964;

    //         for (auto &n:N){
    //             timeOfOrderParameterDistribution[n]=std::vector<double>{0.94,0.945,0.95,0.955,0.96,0.965,0.97,0.975,0.98};
    //             orderParameterOfClusterSizeDistribution[n]=std::vector<double>{0.855,0.86,0.865,0.87,0.975,0.88,0.885,0.89,0.895,0.9,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945};
    //         }
    //         // orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.95,0.99};
    //     }
    //     else if (t_g==0.6){
    //         m_c=0.784;
    //         t_c = 0.901;
    //         for (auto &n:N){
    //             timeOfOrderParameterDistribution[n]=std::vector<double>{0.88,0.885,0.89,0.895,0.9,0.905,0.91,0.915,0.92};
    //             orderParameterOfClusterSizeDistribution[n]=std::vector<double>{0.74,0.745,0.755,0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.805,0.81};
    //         }
    //         // orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85};
    //     }
    //     else if (t_g==0.7){
    //         m_c=0.693;
    //         t_c = 0.856;

    //         for (auto &n:N){
    //             timeOfOrderParameterDistribution[n]=std::vector<double>{0.835,0.84,0.845,0.85,0.855,0.86,0.865,0.87,0.875};
    //             orderParameterOfClusterSizeDistribution[n]=std::vector<double>{0.605,0.61,0.615,0.62,0.625,0.63,0.635,0.64,0.645,0.655,0.66,0.665,0.67,0.675,0.68,0.685,0.69,0.695};
    //         }
    //         // orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85};
    //     }
    //     else if (t_g==0.9){
    //         m_c=0.415;
    //         t_c = 0.718;

    //         for (auto &n:N){
    //             timeOfOrderParameterDistribution[n]=std::vector<double>{0.695,0.7,0.705,0.71,0.715,0.72,0.725,0.73,0.735};
    //             orderParameterOfClusterSizeDistribution[n]=std::vector<double>{0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49};
    //         }
    //         // orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85};
    //     }
    //     else if (t_g==1){
    //         m_c=0;
    //         t_c = 0.5;

    //         orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85};
    //     }
    //     else if (t_g==0.8){
    //         m_c=0.58;
    //         t_c = 0.799;

    //         timeOfOrderParameterDistribution[10000]=std::vector<double>{0.85,0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91,0.915,0.92,0.925};
    //         timeOfOrderParameterDistribution[20000]=std::vector<double>{0.835,0.84,0.845,0.85,0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.90,0.905,0.91};
    //         timeOfOrderParameterDistribution[40000]=std::vector<double>{0.82,0.825,0.83,0.835,0.84,0.845,0.85,0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895};
    //         timeOfOrderParameterDistribution[80000]=std::vector<double>{0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85,0.855,0.86,0.865,0.87,0.875,0.88};
    //         timeOfOrderParameterDistribution[160000]=std::vector<double>{0.79,0.795,0.80,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85,0.855,0.86,0.865};
    //         timeOfOrderParameterDistribution[320000]=std::vector<double>{0.79,0.795,0.80,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85,0.855,0.86,0.865};
    //         timeOfOrderParameterDistribution[640000]=std::vector<double>{0.79,0.795,0.80,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85,0.855,0.86,0.865};
    //         timeOfOrderParameterDistribution[1280000]=std::vector<double>{0.77,0.775,0.78,0.785,0.79,0.795,0.80,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845};
    //         timeOfOrderParameterDistribution[2560000]=std::vector<double>{0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.80,0.805,0.81,0.815,0.82};
    //         timeOfOrderParameterDistribution[5120000]=std::vector<double>{0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.80,0.805,0.81,0.815,0.82};
    //         timeOfOrderParameterDistribution[10240000]=std::vector<double>{0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.80,0.805,0.81,0.815,0.82};


    //         orderParameterOfClusterSizeDistribution[10000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[20000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[40000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[80000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[160000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[320000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[640000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         // orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9};
    //         orderParameterOfClusterSizeDistribution[1280000]=std::vector<double>{0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59};
    //         orderParameterOfClusterSizeDistribution[2560000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[5120000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};
    //         orderParameterOfClusterSizeDistribution[10240000]=std::vector<double>{0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7};

    //     }

    //     return std::make_tuple(timeOfOrderParameterDistribution[t_networkSize], orderParameterOfClusterSizeDistribution[t_networkSize], m_c, t_c);
    // }

} //* End of namespace mBFW::parameters

