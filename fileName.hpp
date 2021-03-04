#pragma once

#include <string>

#include "../library/stringFormat.hpp"

namespace mBFW {
//* File name conventions for storing data of mBFW project
namespace fileName {
inline const std::string base(const int& t_networkSize, const double& t_acceptanceThreshold);

const std::string NGE(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_ensembleSize, const int& t_coreNum = -1);

const std::string NGET(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_ensembleSize, const double& t_time, const int& t_coreNum = -1);

const std::string NGEOP(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_ensembleSize, const double& t_orderParameter, const int& t_coreNum = -1);
}  // namespace fileName
}  // namespace mBFW

inline const std::string mBFW::fileName::base(const int& t_networkSize, const double& t_acceptanceThreshold) {
    return "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G" + to_stringWithPrecision(t_acceptanceThreshold, 1);
}

const std::string mBFW::fileName::NGE(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_ensembleSize, const int& t_coreNum) {
    const std::string fileName = base(t_networkSize, t_acceptanceThreshold) + ",E" + std::to_string(t_ensembleSize);
    return t_coreNum == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(t_coreNum) + ".txt";
}

const std::string mBFW::fileName::NGET(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_ensembleSize, const double& t_time, const int& t_coreNum) {
    const std::string fileName = base(t_networkSize, t_acceptanceThreshold) + ",E" + std::to_string(t_ensembleSize) + ",T" + to_stringWithPrecision(t_time, 4);
    return t_coreNum == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(t_coreNum) + ".txt";
}


const std::string mBFW::fileName::NGEOP(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_ensembleSize, const double& t_orderParameter, const int& t_coreNum) {
    const std::string fileName = base(t_networkSize, t_acceptanceThreshold) + ",E" + std::to_string(t_ensembleSize) + ",OP" + to_stringWithPrecision(t_orderParameter, 4);
    return t_coreNum == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(t_coreNum) + ".txt";
}