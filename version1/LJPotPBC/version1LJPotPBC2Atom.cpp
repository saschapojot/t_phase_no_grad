//
// Created by polya on 5/23/24.
//
#include "version1LJPotPBC2Atom.hpp"


void version1dLJPot2Atom::parseCSV(const int &rowNum, double &alpha1, double &beta1, double &p1, double &q1,
                                   double &alpha2, double &beta2, double &p2, double &q2,double &r0) {

    std::string filePath = "./version1Input/1d/LJPotPBC/V1LJ2Atom1d.csv";
    std::string pyFile = "./version1/LJPotPBC/readCSV.py";
    std::string commandToReadCSV = "python3 " + pyFile + " " + filePath + " " + std::to_string(rowNum);

    std::string result = execPython(commandToReadCSV.c_str());

    std::regex pattern_alpha1("alpha1([+-]?\\d+(\\.\\d+)?)beta1");
    std::smatch match_alpha1;
    if (std::regex_search(result, match_alpha1, pattern_alpha1)) {
        alpha1 = std::stod(match_alpha1[1].str());
    }

    std::regex pattern_beta1("beta1([+-]?\\d+(\\.\\d+)?)p1");
    std::smatch match_beta1;
    if (std::regex_search(result, match_beta1, pattern_beta1)) {
        beta1 = std::stod(match_beta1[1].str());
    }

    std::regex pattern_p1("p1([+-]?\\d+(\\.\\d+)?)q1");
    std::smatch match_p1;
    if (std::regex_search(result, match_p1, pattern_p1)) {
        p1 = std::stod(match_p1[1].str());
    }

    std::regex pattern_q1("q1([+-]?\\d+(\\.\\d+)?)alpha2");
    std::smatch match_q1;
    if (std::regex_search(result, match_q1, pattern_q1)) {
        q1 = std::stod(match_q1[1].str());
    }

    std::regex pattern_alpha2("alpha2([+-]?\\d+(\\.\\d+)?)beta2");
    std::smatch match_alpha2;
    if (std::regex_search(result, match_alpha2, pattern_alpha2)) {
        alpha2 = std::stod(match_alpha2[1].str());
    }

    std::regex pattern_beta2("beta2([+-]?\\d+(\\.\\d+)?)p2");
    std::smatch match_beta2;
    if (std::regex_search(result, match_beta2, pattern_beta2)) {
        beta2 = std::stod(match_beta2[1].str());
    }

    std::regex pattern_p2("p2([+-]?\\d+(\\.\\d+)?)q2");
    std::smatch match_p2;
    if (std::regex_search(result, match_p2, pattern_p2)) {
        p2 = std::stod(match_p2[1].str());
    }

    std::regex pattern_q2("q2([+-]?\\d+(\\.\\d+)?)");
    std::smatch match_q2;
    if (std::regex_search(result, match_q2, pattern_q2)) {
        q2 = std::stod(match_q2[1].str());
    }

    std::regex pattern_x("x=([+-]?\\d+(\\.\\d+)?([eE][-+]?\\d+)?)");
    std::smatch match_x;
    if(std::regex_search(result,match_x,pattern_x)){
        r0=std::stod(match_x[1].str());

    }
    if (r0<0){
        r0=1;
    }
//std::cout<<"r0="<<r0<<std::endl;



}

///
/// @param cmd python execution string
/// @return signal from the python
std::string version1dLJPot2Atom::execPython(const char *cmd) {
    std::array<char, 4096*4> buffer; // Buffer to store command output
    std::string result; // String to accumulate output

    // Open a pipe to read the output of the executed command
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    // Read the output a chunk at a time and append it to the result string
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    return result; // Return the accumulated output



}


///
/// @param xA positions of atom A
/// @param xB positions of atom B
/// @param L total length
/// @return beta*potential
double version1dLJPot2Atom::f(const arma::dcolvec &xA, const arma::dcolvec &xB, const double & L){
    return this->beta * ((*potFuncPtr)(xA, xB,L));



}


///
/// @param xACurr positions of atom A
/// @param xBCurr positions of atom B
/// @param LCurr total length
/// @param zANext proposed positions of atom A
/// @param zBNext proposed positions of atom B
/// @param LNext proposed value of length
void version1dLJPot2Atom::proposal(const arma::dcolvec &xACurr, const arma::dcolvec &xBCurr, const double &LCurr,
              arma::dcolvec &zANext, arma::dcolvec &zBNext, double &LNext){

    std::random_device rd;
    std::ranlux24_base gen(rd());
    //fix left end (0A)
//    zANext(0)=xACurr(0);

    for (int j = 0; j < N; j++) {
        std::normal_distribution<double> dTmp(xACurr(j), stddev);
        zANext(j) = dTmp(gen);
    }


    for (int j = 0; j < N; j++) {
        std::normal_distribution<double> dTmp(xBCurr(j), stddev);
        zBNext(j) = dTmp(gen);
    }

    std::normal_distribution<double> dLast(LCurr,stddev);
    LNext=dLast(gen);



}


///
/// @param xA current positions of atom A
/// @param xB current positions of atom B
/// @param LCurr total length
/// @param zA proposed positions of atom A
/// @param zB proposed positions of atom B
/// @param LNext proposed value of length
/// @return
double version1dLJPot2Atom::acceptanceRatio(const arma::dcolvec &xA, const arma::dcolvec &xB,const double& LCurr,
                       const arma::dcolvec &zA, const arma::dcolvec &zB, const double &LNext){


    double numerator = -f(zA, zB, LNext);

    double denominator = -f(xA, xB,LCurr);
//    double UCurr=(*potFuncPtr)(xA, xB,LCurr);
//    double UNext=(*potFuncPtr)(zA, zB, LNext);
//    std::cout<<"UCurr="<<UCurr<<", UNext="<<UNext<<std::endl;
    double ratio = std::exp(numerator - denominator);

    return std::min(1.0, ratio);

}



///
/// @param xAInit initial positions of A
/// @param xBInit initial positions of B
/// @param LInit
void version1dLJPot2Atom::initPositionsEquiDistance(arma::dcolvec &xAInit, arma::dcolvec &xBInit, double &LInit){
    double a = 5;
    for (int n = 0; n < N; n++) {
        double nDB = static_cast<double >(n);
        xAInit(n) = a * nDB * 2.0;
        xBInit(n) = a * (2.0 * nDB + 1);
    }

    LInit=xBInit(N-1)+2*a;



}


///
/// @param lag decorrelation length
/// @param loopTotal total mc steps
/// @param equilibrium whether equilibrium has reached
/// @param same whether all values of potential are the same
/// @param xALast last positions of atom A
/// @param xBLast last positions of atom B
/// @param LLast last value of total length
void version1dLJPot2Atom::readEqMc(int &lag, int &loopTotal, bool &equilibrium, bool &same, std::vector<double> &xALast,
              std::vector<double> &xBLast, double &LLast){

    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)

    arma::dcolvec xACurr(N);
    arma::dcolvec xBCurr(N);
    double LCurr=0;

    this->initPositionsEquiDistance(xACurr, xBCurr,LCurr);
    double UCurr = (*potFuncPtr)(xACurr, xBCurr,LCurr);

    //output directory
    std::ostringstream sObjT;
    sObjT << std::fixed;
    sObjT << std::setprecision(10);
    sObjT << T;
    std::string TStr = sObjT.str();
    std::string funcName = demangle(typeid(*potFuncPtr).name());
//    std::string  initFuncName= demangle(typeid(initFuncName).name());
    std::string outDir = "./version1Data/1d/func" + funcName +"/row"+std::to_string(rowNum)+ "/T" + TStr + "/";



}