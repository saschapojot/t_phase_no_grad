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
//this subroutine computes the mc evolution for a 1d system, 2-atom, quadratic potential +PBC

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
    virtual double operator()(const double&r,const double &theta0A, const double &theta0B, const double& theta1A, const double& theta1B) const = 0;
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


class quadraticInsideL : public potentialFunction {
public:
    quadraticInsideL(const double &a1Val, const double &a2Val, const double& coef1Val, const double& coef2Val, const double &mAVal, const double &mBVal):potentialFunction(a1Val, a2Val, coef1Val, coef2Val, mAVal, mBVal)  {


        this->a1=a1Val;
        this->a2=a2Val;

        this->coef1=coef1Val;
        this->coef2=coef2Val;

        this->mA=mAVal;
        this->mB=mBVal;

//        std::cout<<"a1="<<this->a1<<", a2="<<this->a2<<", c1="<<this->coef1<<", c2="<<coef2<<std::endl;

    }//end of constructor

public:

    double operator()(const double&r,const double &theta0A, const double &theta0B, const double& theta1A, const double& theta1B) const override {
        double d0A0B=r*(theta0B-theta0A);
        double d0B1A=r*(theta1A-theta0B);
        double d1A1B=r*(theta1B-theta1A);

        double d1B0A=r*(2*PI+theta0A-theta1B);

        double val=coef1*std::pow(d0A0B-a1,2)+coef2*std::pow(d0B1A-a2,2)
                   +coef1*std::pow(d1A1B-a1,2)+coef2*std::pow(d1B0A-a2,2);
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



class version1InsideLQuadratic {
public:
    version1InsideLQuadratic (int rowNum, double temperature, unsigned long long cellNum,
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
        this->h=0.5;//stepSize;


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
    /// @param r
    /// @param theta0A
    /// @param theta0B
    /// @param theta1A
    /// @param theta1B
    /// @return
    double f(const double &r,const double& theta0A, const double &theta0B, const double&theta1A, const double&theta1B);

    ///
    /// @param rCurr
    /// @param sigma
    /// @return
    static double generate_nearby_positive_value(const double& rCurr, const double& sigma) {
        std::random_device rd;  // Random number generator
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::normal_distribution<> d(rCurr, sigma); // Normal distribution with mean rCurr and standard deviation sigma

        double rNext;
        do {
            rNext = d(gen);
        } while (rNext <= 0); // Ensure the generated value is positive

        return rNext;
    }

    static double generate_nearby_0_2pi(const double& theta, const double& sigma){
        std::random_device rd;  // Random number generator
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::normal_distribution<> d(theta, sigma); // Normal distribution with mean rCurr and standard deviation sigma
        double thetaNext;
        do {
            thetaNext = d(gen);
//            std::cout<<"generated value/2pi: "<<thetaNext<<std::endl;
        } while (thetaNext < 0 or thetaNext>2*PI); // Ensure the generated value is positive

        return thetaNext;

    }
    ///
    /// @param rCurr
    /// @param theta0ACurr
    /// @param theta0BCurr
    /// @param theta1ACurr
    /// @param theta1BCurr
    /// @param rNext
    /// @param theta0ANext
    /// @param theta0BNext
    /// @param theta1ANext
    /// @param theta1BNext
    void proposal(const double &rCurr, const double& theta0ACurr,const double&theta0BCurr, const double& theta1ACurr,const double &theta1BCurr,
                  double & rNext, double &theta0ANext, double &theta0BNext, double &theta1ANext, double &theta1BNext);


    double acceptanceRatio(const double &rCurr,const double &theta0ACurr, const double &theta0BCurr, const double& theta1ACurr, const double &theta1BCurr,
                           const double &rNext, const double&theta0ANext, const double &theta0BNext, const double &theta1ANext, const double &theta1BNext);


    ///
    /// @param rInit
    /// @param theta0AInit
    /// @param theta0BInit
    /// @param theta1AInit
    /// @param theta1BInit
    void initPositionsEquiDistance(double &rInit,double &theta0AInit, double &theta0BInit,  double &theta1AInit, double &theta1BInit);

    ///
    /// @param lag
    /// @param loopEq
    /// @param equilibrium
    /// @param same
    /// @param rLast
    /// @param theta0ALast
    /// @param theta0BLast
    /// @param theta1ALast
    /// @param theta1BLast
    /// @param U_ptr
    /// @param r_ptr
    /// @param theta0A_ptr
    /// @param theta0B_ptr
    /// @param theta1A_ptr
    /// @param theta1B_ptr
    void readEqMc(unsigned long long &lag,  unsigned long long &loopTotal,bool &equilibrium, bool& same, double &rLast,
                 double &theta0ALast, double &theta0BLast,  double &theta1ALast, double &theta1BLast,double * U_ptr,double *r_ptr, double *theta0A_ptr,double *theta0B_ptr, double *theta1A_ptr, double * theta1B_ptr);




    ///
    /// @param lag
    /// @param loopEq
    /// @param rInit
    /// @param theta0AInit
    /// @param theta0BInit
    /// @param theta1AInit
    /// @param theta1BInit
    /// @param U_ptr
    /// @param r_ptr
    /// @param theta0A_ptr
    /// @param theta0B_ptr
    /// @param theta1A_ptr
    /// @param theta1B_ptr
    void executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const double &rInit,
                           const double &theta0AInit, const double &theta0BInit, const double &theta1AInit, const double &theta1BInit,double * U_ptr,double *r_ptr,double *theta0A_ptr, double *theta0B_ptr, double *theta1A_ptr, double * theta1B_ptr);


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