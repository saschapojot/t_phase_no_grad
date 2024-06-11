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
    virtual double operator()(const double&L,const double &x0A, const double &x0B, const double& x1A, const double& x1B) const = 0;
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


class quadraticCartesian : public potentialFunction {
public:
    quadraticCartesian(const double &a1Val, const double &a2Val, const double& coef1Val, const double& coef2Val, const double &mAVal, const double &mBVal):potentialFunction(a1Val, a2Val, coef1Val, coef2Val, mAVal, mBVal)  {


        this->a1=a1Val;
        this->a2=a2Val;

        this->coef1=coef1Val;
        this->coef2=coef2Val;

        this->mA=mAVal;
        this->mB=mBVal;

//        std::cout<<"a1="<<this->a1<<", a2="<<this->a2<<", c1="<<this->coef1<<", c2="<<coef2<<std::endl;

    }//end of constructor

public:

    double operator()(const double&L,const double &x0A, const double &x0B, const double& x1A, const double& x1B) const override {
        double d0A0B=x0B-x0A ;
        double d0B1A=x1A-x0B ;
        double d1A1B=x1B-x1A;

        double d1B0A=x0A-x1B+L;

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



class version1CartesianQuadratic {
public:
    version1CartesianQuadratic (int rowNum, double temperature, unsigned long long cellNum,
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
    /// @param x0A
    /// @param x0B
    /// @param x1A
    /// @param x1B
    /// @return
    double f(const double &L,const double& x0A, const double &x0B, const double&x1A, const double&x1B);

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

    static double generate_nearby_0_L(const double& x, const double& sigma,const double &L){
        std::random_device rd;  // Random number generator
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::normal_distribution<> d(x, sigma); // Normal distribution with mean rCurr and standard deviation sigma
        double xNext;
        do {
            xNext = d(gen);
        } while (xNext < 0 or xNext>L); // Ensure the generated value is positive

        return xNext;

    }

    ///
    /// @param LCurr
    /// @param x0ACurr
    /// @param x0BCurr
    /// @param x1ACurr
    /// @param x1BCurr
    /// @param LNext
    /// @param x0ANext
    /// @param x0BNext
    /// @param x1ANext
    /// @param x1BNext
    void proposal(const double &LCurr, const double& x0ACurr,const double&x0BCurr, const double& x1ACurr,const double &x1BCurr,
                  double & LNext, double &x0ANext, double &x0BNext, double &x1ANext, double &x1BNext);

    ///
    /// @param LCurr
    /// @param x0ACurr
    /// @param x0BCurr
    /// @param x1ACurr
    /// @param x1BCurr
    /// @param LNext
    /// @param x0ANext
    /// @param x0BNext
    /// @param x1ANext
    /// @param x1BNext
    /// @return
    double acceptanceRatio(const double &LCurr,const double &x0ACurr, const double &x0BCurr, const double& x1ACurr, const double &x1BCurr,
                           const double &LNext, const double&x0ANext, const double &x0BNext, const double &x1ANext, const double &x1BNext);


   ///
   /// @param LInit
   /// @param x0AInit
   /// @param x0BInit
   /// @param x1AInit
   /// @param x1BInit
    void initPositionsEquiDistance(double &LInit,double &x0AInit, double &x0BInit,  double &x1AInit, double &x1BInit);

    ///
    /// @param lag
    /// @param loopTotal
    /// @param equilibrium
    /// @param same
    /// @param LLast
    /// @param x0ALast
    /// @param x0BLast
    /// @param x1ALast
    /// @param x1BLast
    /// @param U_ptr
    /// @param L_ptr
    /// @param x0A_ptr
    /// @param x0B_ptr
    /// @param x1A_ptr
    /// @param x1B_ptr
    void readEqMc(unsigned long long &lag,  unsigned long long &loopTotal,bool &equilibrium, bool& same, double &LLast,
                 double &x0ALast, double &x0BLast,  double &x1ALast, double &x1BLast,double * U_ptr,double *L_ptr, double *x0A_ptr,double *x0B_ptr, double *x1A_ptr, double * x1B_ptr);




    ///
    /// @param lag
    /// @param loopEq
    /// @param LInit
    /// @param x0AInit
    /// @param x0BInit
    /// @param x1AInit
    /// @param x1BInit
    /// @param U_ptr
    /// @param L_ptr
    /// @param x0A_ptr
    /// @param x0B_ptr
    /// @param x1A_ptr
    /// @param x1B_ptr
    void executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const double &LInit,
                           const double &x0AInit, const double &x0BInit, const double &x1AInit, const double &x1BInit,double * U_ptr,double *L_ptr,double *x0A_ptr, double *x0B_ptr, double *x1A_ptr, double * x1B_ptr);


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