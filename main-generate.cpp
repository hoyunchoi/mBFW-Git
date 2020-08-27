#include <chrono>

#include "generate.hpp"

int main(int argc, char *argv[]){
    //* Input parameters
    const int networkSize=std::stoul(argv[1]);
    const double g=std::stod(argv[2]);
    const int ensembleSize=std::stoul(argv[3]);
    const std::string machine=argv[4];
    const int coreNum=std::stoul(argv[5]);
    constexpr int randomEngineSeed = -1;

    //* Set precision
    double precision;
    g==0.2 ? precision=1e3 : precision=1e4;

    //* Determine which observables to calculate
    std::vector<bool> observables(15);
    observables[0] = true;      //! Order Parameter
    observables[1] = true;      //! Mean Cluster Size
    observables[2] = true;      //! Second Giant
    observables[3] = true;      //! Inter Event Time
    observables[4] = true;      //! Delta Acceptance
    observables[5] = true;      //! Order Parameter Distribution
    observables[6] = true;      //! Cluster Size Distribution
    observables[7] = true;      //! Age Distribution
    observables[8] = true;      //! Inter Event Time Distribution
    observables[9] = true;      //! Delta Upper Bound Distribution
    observables[10] = true;     //! Delta Acceptance Distribution
    observables[11] = true;     //! Inter Event Time vs Delta Acceptance
    observables[12] = true;     //! Upper Bound vs Delta Acceptance
    observables[13] = true;     //! Delta Upper Bound vs Delta Acceptance
    observables[14] = false;    //! Dynamics

    //* run mBFW
    auto start=std::chrono::system_clock::now();
    mBFW::generate::setParameters(networkSize, ensembleSize, g, precision, coreNum, randomEngineSeed, observables);
    mBFW::generate::run();
    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf(" %.6fs for N=%.1e, g=%.1f, ensemble=%d-%d at %s\n", sec.count(),(double)networkSize, g, ensembleSize, coreNum, machine.c_str());

    //* save parameters
    start = std::chrono::system_clock::now();
    mBFW::generate::save();
    sec=std::chrono::system_clock::now()-start;
    printf(" %0.6fs for saving\n", sec.count());

    return 0;
}