//
// Created by polya on 5/23/24.
//

#ifndef T_PHASE_NO_GRAD_VERSION1LJPOTPBC2ATOM_HPP
#define T_PHASE_NO_GRAD_VERSION1LJPOTPBC2ATOM_HPP
#include <algorithm>
#include <armadillo>
#include <array>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/json.hpp>
#include <boost/python.hpp>
#include <boost/python/object/pickle_support.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>


#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cxxabi.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <msgpack.hpp>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

namespace fs = boost::filesystem;

//this subroutine computes the mc evolution for a 1d system, 2-atom, Lennard-Jones+quartic potential +PBC

class potentialFunction {
    //base class for potential function
public:
    potentialFunction (const double &alpha1Val, const double &alpha2Val, const double &beta1Val,
             const double &beta2Val, const double &p1Val, const double &p2Val, const double &q1Val, const double &q2Val, const double &r0) {

        this->alpha1 = alpha1Val;
        this->alpha2 = alpha2Val;
        this->beta1 = beta1Val;
        this->beta2 = beta2Val;
        this->p1 = p1Val;
        this->p2 = p2Val;
        this->q1 = q1Val;
        this->q2 = q2Val;
        this->r0=r0;

    }//end of constructor
    virtual double operator()(const arma::dcolvec &xA, const arma::dcolvec &xB, const double &L) const = 0;
    virtual double dVEst(const double &r, const int &N)const = 0;
    virtual ~ potentialFunction() {};

public:
    double alpha1 = 0;
    double alpha2 = 0;
    double beta1 = 0;
    double beta2 = 0;
    double p1 = 0;
    double p2 = 0;
    double q1 = 0;
    double q2 = 0;
    double r0=0;// eq distance

};


class LJPotPBC : public potentialFunction {
public:
    LJPotPBC(const double &alpha1Val, const double &alpha2Val, const double &beta1Val,
             const double &beta2Val, const double &p1Val, const double &p2Val, const double &q1Val, const double &q2Val, const double &r0Val):potentialFunction(alpha1Val, alpha2Val, beta1Val, beta2Val, p1Val, p2Val, q1Val, q2Val, r0Val)  {}//end of constructor

public:
    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @return potential energy
    double operator()(const arma::dcolvec &xA, const arma::dcolvec &xB, const double& L) const override {
        return V1Total(xA, xB) + V2Total(xA, xB,L);

    }


    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @return the sum of all V1 energy
    double V1Total(const arma::dcolvec &xA, const arma::dcolvec &xB) const {
        arma::dcolvec rVec = arma::abs(xA - xB);


        arma::dcolvec vecPart1 = alpha1 * arma::pow(rVec, -p1);
        arma::dcolvec vecPart2 = -beta1 * arma::pow(rVec, -q1);
        arma::dcolvec vecPart3 = arma::pow(rVec, 4);

        double val = arma::sum(vecPart1) + arma::sum(vecPart2) + arma::sum(vecPart3);
//        std::cout<<"V1Total="<<val<<std::endl;
        return val;


    }

    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @param L length of the PBC loop
    /// @return the sum of all V2 energy under OBC
    double V2Total(const arma::dcolvec &xA, const arma::dcolvec &xB, const double &L) const {
        int N = xB.size();
        if (N <= 1) {
            return 0;
        }
        arma::dcolvec sliceA = xA.subvec(1, N - 1);
//    std::cout<<"sliceA="<<sliceA<<std::endl;
        arma::dcolvec sliceB = xB.subvec(0, N - 2);
        arma::dcolvec rVec = arma::abs(sliceA - sliceB);
//        std::cout<<"sliceB="<<sliceB<<std::endl;

//        std::cout<<"V2Total: rVec="<<rVec<<std::endl;

        arma::dcolvec vecPart1 = alpha2 * arma::pow(rVec, -p2);
        arma::dcolvec vecPart2 = -beta2 * arma::pow(rVec, -q2);
        arma::dcolvec vecPart3 = arma::pow(rVec, 4);

        double rLastAbs=std::abs(L-(xB(N-1)-xA(0)));
        double valBoundary=alpha2*std::pow(rLastAbs,-p2)-beta2*std::pow(rLastAbs,-q2)+std::pow(rLastAbs,4);

        double val = arma::sum(vecPart1) + arma::sum(vecPart2) + arma::sum(vecPart3)+valBoundary;
//        std::cout<<"V2Total="<<val<<std::endl;

        return val;

    }

    double dV1(const double &r)const{
        return -alpha1*p1*std::pow(r,(-p1-1))+beta1*q1*std::pow(r,(-q1-1))
        +4*std::pow(r,3);

    }

    double dV2(const double &r)const{


        return -alpha2*p2*std::pow(r,(-p2-1))+beta2*q2*std::pow(r,(-q2-1))
        +4*std::pow(r,3);
    }

    double dVEst(const double &r, const int &N)const{
        double val=static_cast<double>(N)*(dV1(r)+ dV2(r));
        return val;

    }

public:
    double alpha1 = 0;
    double alpha2 = 0;
    double beta1 = 0;
    double beta2 = 0;
    double p1 = 0;
    double p2 = 0;
    double q1 = 0;
    double q2 = 0;
    double r0=0;// eq distance

};


class version1dLJPot2Atom {
public:
    version1dLJPot2Atom(int rowNum, double temperature, unsigned long long cellNum,
                        const std::shared_ptr<potentialFunction> &funcPtr) {

        this->rowNum = rowNum;
        this->T = temperature;
        this->beta = 1 / T;
        this->potFuncPtr = funcPtr;

//        this->diag=isDiag;
        this->N = cellNum;

        //estimate step size

        double rEst=funcPtr->r0*2.0;

//        std::cout<<"rEst="<<rEst<<std::endl;
        double dValEst=0.1;

        double stepSize=dValEst*T/(std::abs(funcPtr->dVEst(rEst,cellNum)));
        this->h=stepSize*0.5;

        std::cout<<"h="<<h<<std::endl;
        this->stddev =h;




    }





public:

    std::string demangle(const char *name) {
        int status = -1;
        char *demangled = abi::__cxa_demangle(name, NULL, NULL, &status);
        std::string result(name);
        if (status == 0) {
            result = demangled;
        }
        std::free(demangled);
        return result;
    }

    /// @param rowNum row number
    static void parseCSV(const int &rowNum, double &alpha1, double &beta1, double &p1, double &q1,
                         double &alpha2, double &beta2, double &p2, double &q2, double &r0);

    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char *cmd);


    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @param L total length
    /// @return beta*potential
    double f(const arma::dcolvec &xA, const arma::dcolvec &xB, const double & L);


    ///
    /// @param xACurr positions of atom A
    /// @param xBCurr positions of atom B
    /// @param LCurr total length
    /// @param zANext proposed positions of atom A
    /// @param zBNext proposed positions of atom B
    /// @param LNext proposed value of length
    void proposal(const arma::dcolvec &xACurr, const arma::dcolvec &xBCurr, const double &LCurr,
                  arma::dcolvec &zANext, arma::dcolvec &zBNext, double &LNext);

    ///
    /// @param xA current positions of atom A
    /// @param xB current positions of atom B
    /// @param LCurr total length
    /// @param zA proposed positions of atom A
    /// @param zB proposed positions of atom B
    /// @param LNext proposed value of length
    /// @return
    double acceptanceRatio(const arma::dcolvec &xA, const arma::dcolvec &xB,const double& LCurr,
                           const arma::dcolvec &zA, const arma::dcolvec &zB, const double &LNext);


    ///
    /// @param xAInit initial positions of A
    /// @param xBInit initial positions of B
    /// @param LInit
    void initPositionsEquiDistance(arma::dcolvec &xAInit, arma::dcolvec &xBInit, double &LInit);

    ///
    /// @param lag decorrelation length
    /// @param loopTotal total mc steps
    /// @param equilibrium whether equilibrium has reached
    /// @param same whether all values of potential are the same
    /// @param xALast last positions of atom A
    /// @param xBLast last positions of atom B
    /// @param LLast last value of total length
    void readEqMc(unsigned long long &lag, unsigned long long &loopTotal, bool &equilibrium, bool &same, arma::dcolvec &xALast,
                  arma::dcolvec &xBLast, double &LLast);




    ///
    /// @param lag decorrelation length
    /// @param loopEq total loop numbers in reaching equilibrium
    /// @param xA_init xA from readEqMc
    /// @param xB_init xB from readEqMc
    /// @param LInit L from readEqMc
    void executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const arma::dcolvec &xA_init,
                            const arma::dcolvec &xB_init, const double &LInit);


public:
    double T;// temperature
    double beta;
//    int moveNumInOneFlush = 3000;// flush the results to python every moveNumInOneFlush iterations
//    int flushMaxNum = 7000;
    unsigned long long loopMax=3000*80000;//max number of loop to reach equilibrium
    unsigned long  long loopToWrite=3000*20000;
    unsigned long long dataNumTotal = 2000;
    unsigned long long dataNumInEq=0;
    double h;// step size
    unsigned long long N;//number of unit cells

//    double lastFileNum = 0;
    std::shared_ptr<potentialFunction> potFuncPtr;
    double stddev;
    int rowNum;
    unsigned long long nEqCounterStart=0;// loop number when equilibrium is reached

    unsigned long long writeInterval=100000000;
};


void save_array_to_pickle(const std::unique_ptr<double[]>& ptr, std::size_t size, const std::string& filename);

///to msgpack bin file
void save_to_bin_file(const std::unique_ptr<double[]>& data, unsigned long long size, const std::string& filename);
#endif //T_PHASE_NO_GRAD_VERSION1LJPOTPBC2ATOM_HPP
