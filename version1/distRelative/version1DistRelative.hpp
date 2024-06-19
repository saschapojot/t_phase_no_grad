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
#include <initializer_list>
#include <iostream>
#include <math.h>
#include <memory>
#include <msgpack.hpp>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

namespace fs = boost::filesystem;
const auto PI=M_PI;
//Using angular coordinates, this subroutine computes the mc evolution for a 1d system, 2-atom, quadratic potential +PBC

class potentialFunction {
    //base class for potential function
public:
    potentialFunction (const double &a1Val, const double &a2Val, const double& coef1Val, const double& coef2Val, const double &mAVal, const double &mBVal) {
        this->a1=a1Val;
        this->a2=a2Val;

        this->coef1=coef1Val;
        this->coef2=coef2Val;

        this->mA=mAVal;
        this->mB=mBVal;

    }//end of constructor
    virtual double operator()(const double&L,const double &y0, const double &z0, const double& y1) const = 0;
    virtual double dVEst(const double &r, const unsigned long long &N)const = 0;
    virtual ~ potentialFunction() {};

public:
    double a1 ;
    double a2 ;



    double coef1;
    double coef2;

    double mA;
    double mB;

};


class quadraticDistRelative : public potentialFunction {
public:
    quadraticDistRelative(const double &a1Val, const double &a2Val, const double& coef1Val, const double& coef2Val, const double &mAVal, const double &mBVal):potentialFunction(a1Val, a2Val, coef1Val, coef2Val, mAVal, mBVal)  {


        this->a1=a1Val;
        this->a2=a2Val;

        this->coef1=coef1Val;
        this->coef2=coef2Val;

        this->mA=mAVal;
        this->mB=mBVal;

//        std::cout<<"a1="<<this->a1<<", a2="<<this->a2<<", c1="<<this->coef1<<", c2="<<coef2<<std::endl;

    }//end of constructor

public:

    double operator()(const double&L,const double &y0, const double &z0, const double& y1) const override {


        double val=coef1*std::pow(y0-a1,2)+coef2*std::pow(z0-a2,2)
                   +coef1*std::pow(y1-a1,2)+coef2*std::pow(-y0-z0-y1+L-a2,2);
//        std::cout<<"val="<<val<<std::endl;
        return val;

    }





    double dVEst(const double &r, const unsigned long long &N)const{
        double val=0;
        return val;

    }

public:
    double a1 ;
    double a2 ;

//    double r0;// eq distance

    double coef1;
    double coef2;
    double mA;
    double mB;

};



class version1DistRelative {
public:
    version1DistRelative (int rowNum, double temperature, unsigned long long cellNum,
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
double stepForT1=0.1;
        this->h= stepForT1*T;//stepSize;


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
    static void parseCSV(const int &rowNum, double &a1, double &a2, double &c1, double &c2,double &mA,double&mB);

    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char *cmd);


    ///
    /// @param L
    /// @param y0
    /// @param z0
    /// @param y1
    /// @return
    double f(const double &L,const double& y0, const double &z0, const double&y1);

    ///
    /// @param LCurr
    /// @param sigma
    /// @return
    static double generate_nearby_positive_value(const double& LCurr, const double& sigma) {
        std::random_device rd;  // Random number generator
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::normal_distribution<> d(LCurr, sigma); // Normal distribution with mean rCurr and standard deviation sigma

        double LNext;
        do {
            LNext = d(gen);
        } while (LNext <= 0); // Ensure the generated value is positive

        return LNext;
    }

    static double generate_nearby_mL_L(const double& x, const double& sigma,const double &L){
        std::random_device rd;  // Random number generator
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::normal_distribution<> d(x, sigma); // Normal distribution with mean rCurr and standard deviation sigma
        double xNext;
        do {
            xNext = d(gen);
        } while (xNext < -L or xNext>L); // Ensure the generated value is in [-L,L]

        return xNext;

    }

    ///
    /// @param LCurr
    /// @param y0Curr
    /// @param z0Curr
    /// @param y1Curr
    /// @param LNext
    /// @param y0Next
    /// @param z0Next
    /// @param y1Next
    void proposal(const double &LCurr, const double& y0Curr,const double& z0Curr, const double& y1Curr,
                  double & LNext, double & y0Next, double & z0Next, double & y1Next);

    ///
    /// @param LCurr
    /// @param y0Curr
    /// @param z0Curr
    /// @param y1Curr
    /// @param LNext
    /// @param y0Next
    /// @param z0Next
    /// @param y1Next
    /// @return
    double acceptanceRatio(const double &LCurr,const double &y0Curr, const double &z0Curr, const double& y1Curr,
                           const double &LNext, const double& y0Next, const double & z0Next, const double & y1Next);


   ///
   /// @param LInit
   /// @param y0Init
   /// @param z0Init
   /// @param y1Init
    void initPositionsEquiDistance(double &LInit,double &y0Init, double &z0Init,  double &y1Init);


    ///
    /// @param lag
    /// @param loopTotal
    /// @param equilibrium
    /// @param same
    /// @param LLast
    /// @param y0Last
    /// @param z0Last
    /// @param y1Last
    /// @param U_ptr
    /// @param L_ptr
    /// @param y0_ptr
    /// @param z0_ptr
    /// @param y1_ptr
    void readEqMc(unsigned long long &lag,  unsigned long long &loopTotal,bool &equilibrium, bool& same, double &LLast,
                 double &y0Last, double & z0Last,  double & y1Last,double * U_ptr,double *L_ptr, double *y0_ptr,double *z0_ptr, double *y1_ptr);




    ///
    /// @param lag
    /// @param loopEq
    /// @param LInit
    /// @param y0Init
    /// @param z0Init
    /// @param y1Init
    /// @param U_ptr
    /// @param L_ptr
    /// @param y0_ptr
    /// @param z0_ptr
    /// @param y1_ptr
    void executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const double &LInit,
                           const double &y0Init, const double &z0Init, const double &y1Init,double * U_ptr,double *L_ptr,double *y0_ptr, double *z0_ptr, double *y1_ptr);


public:
    double T;// temperature
    double beta;
//    int moveNumInOneFlush = 3000;// flush the results to python every moveNumInOneFlush iterations
//    int flushMaxNum = 7000;

    static const unsigned long long loopMax=400000000;//max number of loop to reach equilibrium
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