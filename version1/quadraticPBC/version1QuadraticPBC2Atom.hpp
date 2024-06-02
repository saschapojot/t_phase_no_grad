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

//this subroutine computes the mc evolution for a 1d system, 2-atom, quadratic potential +PBC

class potentialFunction {
    //base class for potential function
public:
    potentialFunction (const double &a1Val, const double &a2Val, const double& coef1Val, const double& coef2Val, const double &r0Val) {
    this->a1=a1Val;
    this->a2=a2Val;
    this->r0=r0Val;
    this->coef1=coef1Val;
    this->coef2=coef2Val;

    }//end of constructor
    virtual double operator()(const arma::dcolvec &xA, const arma::dcolvec &xB, const double &L) const = 0;
    virtual double dVEst(const double &r, const unsigned long long &N)const = 0;
    virtual ~ potentialFunction() {};

public:
    double a1 ;
    double a2 ;

    double r0;// eq distance

    double coef1;
    double coef2;

};


class quadratic : public potentialFunction {
public:
    quadratic(const double &a1Val, const double &a2Val, const double& coef1Val, const double& coef2Val, const double &r0Val):potentialFunction(a1Val, a2Val,  coef1Val,  coef2Val, r0Val)  {


        this->a1=a1Val;
        this->a2=a2Val;
        this->r0=r0Val;
        this->coef1=coef1Val;
        this->coef2=coef2Val;

//        std::cout<<"a1="<<this->a1<<", a2="<<this->a2<<", c1="<<this->coef1<<", c2="<<coef2<<std::endl;

    }//end of constructor

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
        arma::dcolvec rVec = xA-xB;
//        rVec.print("rVec");
        rVec+=a1;

//        rVec.print("rVec");


        arma::dcolvec vecPart1 = arma::pow(rVec,2)*coef1;//+arma::pow(rVec,4) ;
//        std::cout<<"vecPart1="<<vecPart1<<std::endl;

        double val = arma::sum(vecPart1);
//        std::cout<<"V1Total="<<val<<std::endl;
        return val;


    }

    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @param L length of the PBC loop
    /// @return the sum of all V2 energy under OBC
    double V2Total(const arma::dcolvec &xA, const arma::dcolvec &xB, const double &L) const {
        size_t N = xB.size();
        if (N == 0) {
            return 0;
        }
        arma::dcolvec sliceA = xA.subvec(1, N - 1);
//    std::cout<<"sliceA="<<sliceA<<std::endl;
        arma::dcolvec sliceB = xB.subvec(0, N - 2);
        arma::dcolvec rVec = sliceB-sliceA;
        rVec+=a2;
//        rVec.print("rVec");
//        std::cout<<"sliceB="<<sliceB<<std::endl;

//        std::cout<<"V2Total: rVec="<<rVec<<std::endl;

        arma::dcolvec vecPart1 = coef2*arma::pow(rVec,2);//+arma::pow(rVec,4);

//        vecPart1.print("vecPart1");
        double V2Last=std::pow(L-(xB(N-1)-xA(0))-a2,2)*coef2;//+std::pow(L-(xB(N-1)-xA(0))-a2,4);

//        std::cout<<"V2Last="<<V2Last<<std::endl;

        double val = arma::sum(vecPart1) + V2Last;
//        std::cout<<"V2Total="<<val<<std::endl;

        return val;

    }




    double dVEst(const double &r, const unsigned long long &N)const{
        double val=0;
        return val;

    }

public:
    double a1 ;
    double a2 ;

    double r0;// eq distance

    double coef1;
    double coef2;

};


class quadraticQuartic : public potentialFunction {
public:
    quadraticQuartic(const double &a1Val, const double &a2Val, const double& coef1Val, const double& coef2Val, const double &r0Val):potentialFunction(a1Val, a2Val,  coef1Val,  coef2Val, r0Val)  {


        this->a1=a1Val;
        this->a2=a2Val;
        this->r0=r0Val;
        this->coef1=coef1Val;
        this->coef2=coef2Val;

//        std::cout<<"a1="<<this->a1<<", a2="<<this->a2<<", c1="<<this->coef1<<", c2="<<coef2<<std::endl;

    }//end of constructor

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
        arma::dcolvec rVec = xA-xB;
//        rVec.print("rVec");
        rVec+=a1;

//        rVec.print("rVec");


        arma::dcolvec vecPart1 = arma::pow(rVec,2)*coef1+20*arma::pow(rVec,4) ;
//        std::cout<<"vecPart1="<<vecPart1<<std::endl;

        double val = arma::sum(vecPart1);
//        std::cout<<"V1Total="<<val<<std::endl;
        return val;


    }

    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @param L length of the PBC loop
    /// @return the sum of all V2 energy under OBC
    double V2Total(const arma::dcolvec &xA, const arma::dcolvec &xB, const double &L) const {
        size_t N = xB.size();
        if (N == 0) {
            return 0;
        }
        arma::dcolvec sliceA = xA.subvec(1, N - 1);
//    std::cout<<"sliceA="<<sliceA<<std::endl;
        arma::dcolvec sliceB = xB.subvec(0, N - 2);
        arma::dcolvec rVec = sliceB-sliceA;
        rVec+=a2;
//        rVec.print("rVec");
//        std::cout<<"sliceB="<<sliceB<<std::endl;

//        std::cout<<"V2Total: rVec="<<rVec<<std::endl;

        arma::dcolvec vecPart1 = coef2*arma::pow(rVec,2)+20*arma::pow(rVec,4);

//        vecPart1.print("vecPart1");
        double V2Last=std::pow(L-(xB(N-1)-xA(0))-a2,2)*coef2+20*std::pow(L-(xB(N-1)-xA(0))-a2,4);

//        std::cout<<"V2Last="<<V2Last<<std::endl;

        double val = arma::sum(vecPart1) + V2Last;
//        std::cout<<"V2Total="<<val<<std::endl;

        return val;

    }




    double dVEst(const double &r, const unsigned long long &N)const{
        double val=0;
        return val;

    }

public:
    double a1 ;
    double a2 ;

    double r0;// eq distance

    double coef1;
    double coef2;

};

class version1Quadratic {
public:
    version1Quadratic (int rowNum, double temperature, unsigned long long cellNum,
                        const std::shared_ptr<potentialFunction> &funcPtr) {

        this->rowNum = rowNum;
        this->T = temperature;
        this->beta = 1 / T;
        this->potFuncPtr = funcPtr;

//        this->diag=isDiag;
        this->N = cellNum;

        //estimate step size


//        double rEst=funcPtr->r0;
//
//        std::cout<<"rEst="<<rEst<<std::endl;
//        double dValEst=2;

//        double stepSize=dValEst*T/(std::abs(funcPtr->dVEst(rEst,cellNum)));
//        if (stepSize>0.005){
//            stepSize=0.005;
//        }
        this->h=0.01;//stepSize;


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
    static void parseCSV(const int &rowNum, double &a1, double &a2, double &c1, double &c2);

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
                  arma::dcolvec &xBLast, double &LLast,double * U_ptr,double *L_ptr, double *xA_ptr, double *xB_ptr );




    ///
    /// @param lag decorrelation length
    /// @param loopEq total loop numbers in reaching equilibrium
    /// @param xA_init xA from readEqMc
    /// @param xB_init xB from readEqMc
    /// @param LInit L from readEqMc
    void executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const arma::dcolvec &xA_init,
                            const arma::dcolvec &xB_init, const double &LInit,double * U_ptr,double *L_ptr, double *xA_ptr, double *xB_ptr );


public:
    double T;// temperature
    double beta;
//    int moveNumInOneFlush = 3000;// flush the results to python every moveNumInOneFlush iterations
//    int flushMaxNum = 7000;

   static const unsigned long long loopMax=100000000;//3000*6000;//max number of loop to reach equilibrium
    static const unsigned long  long loopToWrite=8000000;

    unsigned long long dataNumTotal = 2000;
    unsigned long long dataNumInEq=0;
    double h;// step size
    unsigned long long N;//number of unit cells

//    double lastFileNum = 0;
    std::shared_ptr<potentialFunction> potFuncPtr;
    double stddev;
    int rowNum;
    unsigned long long nEqCounterStart=0;// loop number when equilibrium is reached

    unsigned long long writeInterval=loopToWrite;
};


void save_array_to_pickle(double* ptr, std::size_t size, const std::string& filename);

///to msgpack bin file
void save_to_bin_file(double* data, unsigned long long size, const std::string& filename);
#endif //T_PHASE_NO_GRAD_VERSION1LJPOTPBC2ATOM_HPP
