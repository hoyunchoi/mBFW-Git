# pragma once

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <tuple>

#include "fileName.hpp"

#include "../library/linearAlgebra.hpp"

namespace mBFW::parameters{
    using namespace linearAlgebra;
    const std::string baseDirectory = "../data/mBFW/";
    const std::set<int> networkSizeList = {10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000, 5120000, 10240000};
    const std::set<double> acceptanceThresholdList = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    //* Set parameter for orderParameter_clusterSizeDist
    std::set<double> set_orderParameter_clusterSizeDist(const int& t_networkSize, const double& t_acceptanceThreshold){
        std::set<double> orderParameter_clusterSizeDist;
        if (t_acceptanceThreshold == 0.2){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.992, 0.994, 0.996, 0.998};
            // if (t_networkSize <= 40000){
            //     orderParameter_clusterSizeDist = {0.951, 0.952, 0.953, 0.954, 0.955, 0.956, 0.957, 0.958, 0.959, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969};
            // }
            // else {
            //     orderParameter_clusterSizeDist = {0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969, 0.971, 0.972, 0.973, 0.974, 0.975, 0.976, 0.978, 0.979};
            // }
        }
        else if (t_acceptanceThreshold == 0.3){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.955, 0.96, 0.965, 0.97};
            // orderParameter_clusterSizeDist = {0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.91, 0.92, 0.93, 0.94};
        }
        else if (t_acceptanceThreshold == 0.4){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.905, 0.91, 0.915, 0.92, 0.95};
            // orderParameter_clusterSizeDist = {0.65, 0.66, 0.67, 0.68, 0.69, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.81, 0.82, 0.83, 0.84, 0.85};
            orderParameter_clusterSizeDist = {0.845, 0.855, 0.86, 0.865, 0.870, 0.875, 0.88, 0.885, 0.89, 0.895};
        }
        else if (t_acceptanceThreshold == 0.5){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9};
            // if (t_networkSize <= 20000){
            //     orderParameter_clusterSizeDist = {0.79, 0.791, 0.792, 0.793, 0.794, 0.795, 0.796, 0.797, 0.798, 0.799, 0.801, 0.802, 0.803, 0.804, 0.805, 0.806, 0.807, 0.808, 0.809};
            // }
            // else if (t_networkSize <= 80000){
            //     orderParameter_clusterSizeDist = {0.801, 0.802, 0.803, 0.804, 0.805, 0.806, 0.807, 0.808, 0.809, 0.811, 0.812, 0.813, 0.814, 0.815, 0.816, 0.817, 0.818, 0.819};
            // }
            // else{
            //     orderParameter_clusterSizeDist = {0.811, 0.812, 0.813, 0.814, 0.815, 0.816, 0.817, 0.818, 0.819, 0.821, 0.822, 0.823, 0.824, 0.825, 0.826, 0.827, 0.828, 0.829};
            // }
            if (t_networkSize <= 10000){
                orderParameter_clusterSizeDist = {0.75, 0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785};
            }
        }
        else if (t_acceptanceThreshold == 0.6){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.9};
            // orderParameter_clusterSizeDist = {0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.675, 0.68, 0.685, 0.69, 0.695, 0.705, 0.71};
            if (t_networkSize >= 20000){
                orderParameter_clusterSizeDist = {0.715, 0.725, 0.73, 0.735, 0.745, 0.750, 0.755};
            }
        }
        else if (t_acceptanceThreshold == 0.7){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.67, 0.69, 0.7, 0.71, 0.73, 0.8, 0.9};
            // orderParameter_clusterSizeDist = {0.55, 0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, 0.595, 0.605, 0.61, 0.615, 0.62, 0.625, 0.63, 0.64};
            if (t_networkSize >= 40000){
                orderParameter_clusterSizeDist = {0.635, 0.645, 0.655, 0.660, 0.665};
            }
        }
        else if (t_acceptanceThreshold == 0.8){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.7, 0.8, 0.9};
            // orderParameter_clusterSizeDist = {0.35, 0.36, 0.365, 0.37, 0.375, 0.38, 0.385, 0.39, 0.395, 0.41, 0.415, 0.42, 0.425, 0.43, 0.435, 0.44, 0.45};
            if (t_networkSize == 20000){
                orderParameter_clusterSizeDist = {0.445, 0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495};
            }
            else if (t_networkSize == 40000){
                orderParameter_clusterSizeDist = {0.455, 0.460, 0.465, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495, 0.505, 0.510, 0.515};
            }
            else if (t_networkSize == 80000 || t_networkSize == 160000){
                orderParameter_clusterSizeDist = {0.505, 0.510, 0.515, 0.525, 0.530, 0.535};
            }
            else if (t_networkSize >= 320000){
                orderParameter_clusterSizeDist = {0.525, 0.530, 0.535, 0.545, 0.550, 0.555};
            }
        }
        else if (t_acceptanceThreshold == 0.9){
            // orderParameter_clusterSizeDist = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.5, 0.6, 0.7, 0.8, 0.9};
            // orderParameter_clusterSizeDist = {0.16, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.25};
            if (t_networkSize == 20000 || t_networkSize == 40000){
                orderParameter_clusterSizeDist = {0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.295, 0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.355};
            }
            else if (t_networkSize == 80000){
                orderParameter_clusterSizeDist = {0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.355, 0.365, 0.370, 0.375};
            }
            else if (t_networkSize >= 160000){
                orderParameter_clusterSizeDist = {0.365, 0.370, 0.375, 0.385, 0.390, 0.395};
            }
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
            time_clusterSizeDist = {0.9971, 0.9972, 0.9973, 0.9974, 0.9975, 0.9976, 0.9977, 0.9978, 0.9979};
        }
        else if (t_acceptanceThreshold == 0.3){
            time_clusterSizeDist = {0.95, 0.96, 0.97, 0.98, 0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.989, 0.99, 0.995};
        }
        else if (t_acceptanceThreshold == 0.4){
            time_clusterSizeDist = {0.92, 0.94, 0.96, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969, 0.97, 0.98, 0.99};
        }
        else if (t_acceptanceThreshold == 0.5){
            time_clusterSizeDist = {0.9, 0.91, 0.92, 0.93, 0.931, 0.932, 0.933, 0.934, 0.935, 0.936, 0.937, 0.938, 0.939, 0.94, 0.95, 0.96};
            time_clusterSizeDist = {0.941, 0.942, 0.943, 0.944, 0.945, 0.946, 0.947, 0.948, 0.949};
        }
        else if (t_acceptanceThreshold == 0.6){
            time_clusterSizeDist = {0.84, 0.86, 0.88, 0.89, 0.896, 0.898, 0.90, 0.902, 0.904, 0.905, 0.906, 0.908, 0.91, 0.92, 0.94};
        }
        else if (t_acceptanceThreshold == 0.7){
            time_clusterSizeDist = {0.8, 0.82, 0.83, 0.84, 0.846, 0.848, 0.85, 0.852, 0.854, 0.856, 0.858, 0.86, 0.87};
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

    //* Set parameter for time_orderParameterDist
    std::set<double> set_time_orderParameterDist(const int& t_networkSize, const double& t_acceptanceThreshold){
        std::set<double> time_orderParameterDist;
        if (t_acceptanceThreshold == 0.2){
            time_orderParameterDist = {0.991, 0.992, 0.993, 0.9935, 0.994, 0.9945, 0.995, 0.9955, 0.996, 0.9965, 0.997, 0.998, 0.999};
        }
        else if (t_acceptanceThreshold == 0.3){
            time_orderParameterDist = {0.98, 0.9805, 0.981, 0.9815, 0.982, 0.9825, 0.983, 0.9835, 0.984};
        }
        else if (t_acceptanceThreshold == 0.4){
            time_orderParameterDist = {0.961, 0.9615, 0.962, 0.9625, 0.963, 0.9635, 0.964, 0.9645, 0.965, 0.9655, 0.966};
        }
        else if (t_acceptanceThreshold == 0.5){
            time_orderParameterDist = {0.933, 0.9335, 0.934, 0.9345, 0.935, 0.9355, 0.936, 0.9365, 0.937};
        }
        else if (t_acceptanceThreshold == 0.6){
            time_orderParameterDist = {0.897, 0.8975, 0.898, 0.8985, 0.899, 0.8995, 0.900, 0.9005, 0.901};
        }
        else if (t_acceptanceThreshold == 0.7){
            time_orderParameterDist = {0.852, 0.8525, 0.853, 0.8535, 0.854, 0.8545, 0.855, 0.8555, 0.856, 0.8565, 0.857};
        }
        else if (t_acceptanceThreshold == 0.8){
            time_orderParameterDist = {0.795, 0.7955, 0.796, 0.7965, 0.797, 0.7975, 0.798, 0.7985, 0.799, 0.7995, 0.800};
        }
        else if (t_acceptanceThreshold == 0.9){
            time_orderParameterDist = {0.715, 0.7155, 0.716, 0.7165, 0.717, 0.7175, 0.718, 0.7185, 0.719};
        }
        else if (t_acceptanceThreshold == 1.0){
            time_orderParameterDist = {};
        }
        else{
            std::cout << "Not defined acceptance rate\n";
        }
        return time_orderParameterDist;
    }

    //* Read t_a,m_a,t_c_var,m_c_csd and set
    const std::tuple<double, double, double, double> set_points(const int& t_networkSize, const double& t_acceptanceThreshold){
        double t_a = 0.0;
        double m_a = 0.0;
        double t_c = 0.0;
        double m_c = 0.0;
        std::ifstream readFile(baseDirectory + "points/" + fileName::base(t_networkSize, t_acceptanceThreshold) + ".txt");
        std::string line;
        while (getline(readFile, line)){
            //* Find line for each points
            if (line.find("t_a") != line.npos){
                t_a = std::stod(line.substr(line.find(": ") + 2));
            }
            else if (line.find("m_a") != line.npos){
                m_a = std::stod(line.substr(line.find(": ") + 2));
            }
            else if (line.find("t_c_var") != line.npos){
                t_c = std::stod(line.substr(line.find(": ") + 2));
            }
            else if (line.find("m_c_csd") != line.npos){
                m_c = std::stod(line.substr(line.find(": ") + 2));
            }
        }
        return std::make_tuple(t_a, m_a, t_c, m_c);
    }
} //* End of namespace mBFW::parameters