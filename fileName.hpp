#pragma once

#include "../library/stringFormat.hpp"

namespace mBFW{
    //* File name conventions for storing data of mBFW project
    namespace fileName{
        inline const std::string base(const int& t_networkSize, const double& t_acceptanceThreshold){
            return "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G" + to_stringWithPrecision(t_acceptanceThreshold, 1);
        }

        const std::string NGE(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const int& t_coreNum=-1){
            const std::string fileName = base(t_networkSize, t_acceptanceThreshold) + ",E"+std::to_string(t_ensembleSize);

            //* name for averaged or binned file which doesn't have specific core num
            if (t_coreNum==-1){
                return fileName + ".txt";
            }

            //* name for regular file
            else{
                return fileName + "-" + std::to_string(t_coreNum) + ".txt";
            }
        }

        const std::string NGET(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const double& t_time, const int& t_coreNum=-1){
            const std::string fileName = base(t_networkSize, t_acceptanceThreshold) + ",E"+std::to_string(t_ensembleSize) + ",T"+to_stringWithPrecision(t_time,4);

            return t_coreNum==-1 ? fileName + ".txt" : fileName + "-"+std::to_string(t_coreNum)+".txt";
        }

        const std::string NGEOP(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const double& t_orderParameter, const int& t_coreNum=-1){
            const std::string fileName = "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G"+to_stringWithPrecision(t_acceptanceThreshold,1) + ",E"+std::to_string(t_ensembleSize) + ",OP"+to_stringWithPrecision(t_orderParameter,4);

            return t_coreNum==-1 ? fileName + ".txt" : fileName + "-"+std::to_string(t_coreNum)+".txt";
        }

        const std::string dynamics_txt(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_randomEngineSeed){
            return base(t_networkSize, t_acceptanceThreshold) + "-" + std::to_string(t_randomEngineSeed) + ".txt";
        }
    } //* End of namespace mBFW::fileName
}