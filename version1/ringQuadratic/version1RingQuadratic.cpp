//
// Created by polya on 5/23/24.
//
#include "version1RingQuadratic.hpp"


void version1RingQuadratic::parseCSV(const int &rowNum, double &a1, double &a2, double &c1, double &c2,double &mA, double &mB) {

    std::string filePath = "./version1Input/1d/ringQuadratic/ringQuadratic.csv";
    std::string pyFile = "./version1/ringQuadratic/readCSV.py";
    std::string commandToReadCSV = "python3 " + pyFile + " " + filePath + " " + std::to_string(rowNum);

    std::string result = execPython(commandToReadCSV.c_str());

    std::regex pattern_a1("a1([+-]?\\d+(\\.\\d+)?)a2");
    std::smatch match_a1;
    if (std::regex_search(result, match_a1, pattern_a1)) {
        a1 = std::stod(match_a1[1].str());
    }

    std::regex pattern_a2("a2([+-]?\\d+(\\.\\d+)?)c1");
    std::smatch match_a2;
    if (std::regex_search(result, match_a2, pattern_a2)) {
        a2= std::stod(match_a2[1].str());
    }

    std::regex pattern_c1("c1([+-]?\\d+(\\.\\d+)?)c2");
    std::smatch match_c1;
    if (std::regex_search(result, match_c1, pattern_c1)) {
        c1= std::stod(match_c1[1].str());
    }

    std::regex pattern_c2("c2([+-]?\\d+(\\.\\d+)?)");
    std::smatch match_c2;
    if (std::regex_search(result, match_c2, pattern_c2)) {
       c2 = std::stod(match_c2[1].str());
    }

    std::regex pattern_mA("mA([+-]?\\d+(\\.\\d+)?)mB");
    std::smatch match_mA;
    if (std::regex_search(result, match_mA, pattern_mA)) {
        mA = std::stod(match_mA[1].str());
    }


    std::regex pattern_mB("mB([+-]?\\d+(\\.\\d+)?)");
    std::smatch match_mB;
    if (std::regex_search(result, match_mB, pattern_mB)) {
        mB = std::stod(match_mB[1].str());
    }










//std::cout<<"r0="<<r0<<std::endl;



}

///
/// @param cmd python execution string
/// @return signal from the python
std::string version1RingQuadratic::execPython(const char *cmd) {
    std::array<char, 4096*10> buffer; // Buffer to store command output
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
double version1RingQuadratic::f(const double &r, const double &theta0B, const double&theta1A, const double&theta1B){
    return this->beta * ((*potFuncPtr)(r, theta0B,theta1A,theta1B));



}


///
/// @param rCurr
/// @param theta0BCurr
/// @param theta1ACurr
/// @param theta1BCurr
/// @param rNext
/// @param theta0BNext
/// @param theta1ANext
/// @param theta1BNext
void version1RingQuadratic::proposal(const double &rCurr,const double&theta0BCurr, const double& theta1ACurr,const double &theta1BCurr,
                                     double & rNext, double &theta0BNext, double &theta1ANext, double &theta1BNext){

//    std::random_device rd;
//    std::ranlux24_base gen(rd());

    //next theta0B
//    std::normal_distribution<double> dTmpTheta0B(theta0BCurr,stddev);
//    theta0BNext=dTmpTheta0B(gen);
//    theta0BNext=std::fmod(theta0BNext,2*PI);
theta0BNext= generate_nearby_0_2pi(theta0BCurr,stddev);

    //next theta1A
//    std::normal_distribution<double> dTmpTheta1A(theta1ACurr,stddev);
//    theta1ANext=dTmpTheta1A(gen);
//    theta1ANext=std::fmod(theta1ANext,2*PI);
    theta1ANext= generate_nearby_0_2pi(theta1ACurr,stddev);

    //next theta1B
//    std::normal_distribution<double> dTmpTheta1B(theta1BCurr,stddev);
//    theta1BNext=dTmpTheta1B(gen);
//    theta1BNext=std::fmod(theta1BNext,2*PI);
theta1BNext= generate_nearby_0_2pi(theta1BNext,stddev);

    //next r

    rNext= generate_nearby_positive_value(rCurr,stddev);





}


///
/// @param rCurr
/// @param theta0BCurr
/// @param theta1ACurr
/// @param theta1BCurr
/// @param rNext
/// @param theta0BNext
/// @param theta1ANext
/// @param theta1BNext
/// @return
double version1RingQuadratic::acceptanceRatio(const double &rCurr, const double &theta0BCurr, const double& theta1ACurr, const double &theta1BCurr,
                                              const double &rNext, const double &theta0BNext, const double &theta1ANext, const double &theta1BNext){


    double numerator = -f(rNext, theta0BNext, theta1ANext, theta1BNext);

    double denominator = -f(rCurr,theta0BCurr,theta1ACurr,theta1BCurr);
//    double UCurr=(*potFuncPtr)(xA, xB,LCurr);
//    double UNext=(*potFuncPtr)(zA, zB, LNext);
//    std::cout<<"UCurr="<<UCurr<<", UNext="<<UNext<<std::endl;
    double ratio = std::exp(numerator - denominator);

    return std::min(1.0, ratio);

}



///
/// @param rInit
/// @param theta0BInit
/// @param theta1AInit
/// @param theta1BInit
void version1RingQuadratic::initPositionsEquiDistance(double &rInit, double &theta0BInit,  double &theta1AInit, double &theta1BInit){
    double a = 1.5;
    rInit=2*a;

    theta0BInit=0.9*PI;

    theta1AInit=1*PI;

    theta1BInit=1.5*PI;



}


///
/// @param lag
/// @param loopEq
/// @param rInit
/// @param theta0BInit
/// @param theta1AInit
/// @param theta1BInit
/// @param U_ptr
/// @param r_ptr
/// @param theta0B_ptr
/// @param theta1A_ptr
/// @param theta1B_ptr
void version1RingQuadratic::readEqMc(unsigned long long &lag, unsigned long long &loopTotal,bool &equilibrium, bool& same, double &rLast,
                                     double &theta0BLast,  double &theta1ALast, double &theta1BLast,double * U_ptr,double *r_ptr, double *theta0B_ptr, double *theta1A_ptr, double * theta1B_ptr ){

    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)


    //output directory
    std::ostringstream sObjT;
    sObjT << std::fixed;
    sObjT << std::setprecision(10);
    sObjT << T;
    std::string TStr = sObjT.str();
    std::string funcName = demangle(typeid(*potFuncPtr).name());
//    std::string  initFuncName= demangle(typeid(initFuncName).name());
    std::string outDir = "./version1Data/1d/func" + funcName +"/row"+std::to_string(rowNum)+ "/T" + TStr + "/";
    std::string outUAllPickleSubDir = outDir + "UAllPickle/";
    std::string outUAllBinSubDir = outDir + "UAllBin/";

    std::string out_r_AllBinSubDir = outDir + "r_AllBin/";
    std::string out_theta0B_AllBinSubDir=outDir+"theta0B_AllBin/";

    std::string out_theta1A_AllBinSubDir=outDir+"theta1A_AllBin/";

    std::string out_theta1B_AllBinSubDir=outDir+"theta1B_AllBin/";



    if (!fs::is_directory(outUAllPickleSubDir) || !fs::exists(outUAllPickleSubDir)) {
        fs::create_directories(outUAllPickleSubDir);
    }
    if (!fs::is_directory(outUAllBinSubDir) || !fs::exists(outUAllBinSubDir)) {
        fs::create_directories(outUAllBinSubDir);
    }
    if (!fs::is_directory(out_r_AllBinSubDir) || !fs::exists(out_r_AllBinSubDir)) {
        fs::create_directories(out_r_AllBinSubDir);
    }
    if (!fs::is_directory(out_theta0B_AllBinSubDir ) || !fs::exists(out_theta0B_AllBinSubDir )) {
        fs::create_directories(out_theta0B_AllBinSubDir );
    }

    if (!fs::is_directory(out_theta1A_AllBinSubDir ) || !fs::exists(out_theta1A_AllBinSubDir )) {
        fs::create_directories(out_theta1A_AllBinSubDir );
    }

    if (!fs::is_directory(out_theta1B_AllBinSubDir ) || !fs::exists(out_theta1B_AllBinSubDir )) {
        fs::create_directories(out_theta1B_AllBinSubDir );
    }
//    std::cout<<"loop total="<<loopMax<<std::endl;
    std::regex stopRegex("stop");
    std::regex wrongRegex("wrong");
    std::regex ErrRegex("Err");
    std::regex lagRegex("lag=\\s*(\\d+)");

    std::regex sameRegex("same");
    std::regex eqRegex("equilibrium");
    std::regex ctStartRegex("nCounterStart=\\s*(\\d+)");
    std::regex dataNumEqRegex("numDataPoints=\\s*(\\d+)");

    std::smatch matchUStop;
    std::smatch matchUWrong;
    std::smatch matchUErr;
    std::smatch matchULag;

    std::smatch matchUSame;
    std::smatch matchUEq;
    std::smatch matchCounterStart;
    std::smatch matchDataNumEq;






    double rCurr;
    double theta0BCurr;
    double theta1ACurr;
    double theta1BCurr;


    this->initPositionsEquiDistance(rCurr,theta0BCurr,theta1ACurr,theta1BCurr);
    double UCurr = (*potFuncPtr)(rCurr,theta0BCurr,theta1ACurr,theta1BCurr);
    int lpNum=0;
    bool active = true;
    const auto tMCStart{std::chrono::steady_clock::now()};



    while(lpNum<this->loopMax and active==true){
        //propose a move
        double rNext;
        double theta0BNext;
        double theta1ANext;
        double theta1BNext;

        proposal(rCurr,theta0BCurr,theta1ACurr,theta1BCurr,rNext,theta0BNext,theta1ANext,theta1BNext);
        double r = acceptanceRatio(rCurr,theta0BCurr,theta1ACurr,theta1BCurr,rNext,theta0BNext,theta1ANext,theta1BNext);
        double u = distUnif01(e2);
//        double UTmp=UCurr;
        if (u <= r) {
//                std::cout<<"UCurr="<<UCurr<<std::endl;
            rCurr=rNext;
            theta0BCurr=theta0BNext;
            theta1ACurr=theta1ANext;
            theta1BCurr=theta1BNext;
            UCurr = (*potFuncPtr)(rCurr,theta0BCurr,theta1ACurr,theta1BCurr);
//                std::cout<<"UNext="<<UCurr<<std::endl;
//            double UDiff=UCurr-UTmp;
//            std::cout<<"UDiff="<<UDiff<<std::endl;
//            std::cout<<"-UDiff/T="<<-UDiff/T<<std::endl;
        }//end of accept-reject

        U_ptr[lpNum]=UCurr;
        r_ptr[lpNum]=rCurr;
        theta0B_ptr[lpNum]=theta0BCurr;
        theta1A_ptr[lpNum]=theta1ACurr;
        theta1B_ptr[lpNum]=theta1BCurr;

        //write to file every loopToWrite loops, and inquire equilibrium
        if ((lpNum+1)%loopToWrite==0 and lpNum>1){
            unsigned long long sizeOfArray=lpNum+1;

//            unsigned long long lpStart=0;

            std::string filenameMiddle="loopStart0ReachEq";
            std::string outUPickleFileName=outUAllPickleSubDir+filenameMiddle+ ".UAll.pkl";
            std::string outUBinFileName=outUAllBinSubDir+filenameMiddle+".UAll.bin";
            save_array_to_pickle(U_ptr,sizeOfArray,outUPickleFileName);
            save_to_bin_file(U_ptr,sizeOfArray,outUBinFileName);


            std::string out_r_BinFileName=out_r_AllBinSubDir+filenameMiddle+".r.bin";
            save_to_bin_file(r_ptr,sizeOfArray,out_r_BinFileName);

            std::string out_theta0B_BinFileName=out_theta0B_AllBinSubDir+filenameMiddle+".theta0B_All.bin";
            save_to_bin_file(theta0B_ptr,sizeOfArray,out_theta0B_BinFileName);

            std::string out_theta1A_BinFileName=out_theta1A_AllBinSubDir+filenameMiddle+".theta1A_All.bin";
            save_to_bin_file(theta1A_ptr,sizeOfArray,out_theta1A_BinFileName);

            std::string out_theta1B_BinFileName=out_theta1B_AllBinSubDir+filenameMiddle+".theta1B_All.bin";
            save_to_bin_file(theta1B_ptr,sizeOfArray,out_theta1B_BinFileName);

            const auto tWriteEnd{std::chrono::steady_clock::now()};

            const std::chrono::duration<double> elapsed_seconds{tWriteEnd - tMCStart};
            std::cout << "loop " << lpNum << std::endl;
            std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;

            //communicate with python to inquire equilibrium

            //inquire equilibrium of U
            std::string commandU = "python3 checkVec.py " + outUPickleFileName;
            std::string resultU;

            try {
                const auto tPyStart{std::chrono::steady_clock::now()};
                resultU = this->execPython(commandU.c_str());
                std::cout << "U message from python: " << resultU << std::endl;
                const auto tPyEnd{std::chrono::steady_clock::now()};

                const std::chrono::duration<double> elapsedpy_secondsAll{tPyEnd - tPyStart};
                std::cout << "py time: " << elapsedpy_secondsAll.count()  << " s" << std::endl;

            }
            catch (const std::exception &e) {
                std::cerr << "Error: " << e.what() << std::endl;
                std::exit(10);
            }
            catch (...) {
                // Handle any other exceptions
                std::cerr << "Error" << std::endl;
                std::exit(11);
            }

            // parse result
            if (std::regex_search(resultU, matchUErr, ErrRegex)) {
                std::cout << "error encountered" << std::endl;
                std::exit(12);
            }

            if (std::regex_search(resultU, matchUWrong, wrongRegex)) {
                std::exit(13);
            }

            if (std::regex_search(resultU, matchUStop, stopRegex)) {
                if (std::regex_search(resultU, matchUSame, sameRegex)) {
                    active = false;
                    same = true;



                }


            }//end of regex search for same


            if (std::regex_search(resultU, matchUEq, eqRegex)) {
                if (std::regex_search(resultU, matchULag, lagRegex)) {

                    std::string lagStrU = matchULag.str(1);
                    unsigned long long lagU = std::stoull(lagStrU);
                    std::cout << "lag=" << lagU << std::endl;
                    lag = lagU;


                   if( std::regex_search(resultU,matchCounterStart,ctStartRegex)){
                       this->nEqCounterStart=std::stoull(matchCounterStart.str(1));

                       std::cout<<"nEqCounterStart="<<nEqCounterStart<<std::endl;

                   }


                    if(std::regex_search(resultU,matchDataNumEq,dataNumEqRegex)){
                        this->dataNumInEq=std::stoull(matchDataNumEq.str(1));

                        std::cout<<"dataNumInEq="<<dataNumInEq<<std::endl;

                    }

                    active = false;

                }


            }//end of regex search for equilibrium






        }//end write and inquire (if)




        lpNum++;

    }//end while

    rLast=rCurr;
    theta0BLast=theta0BCurr;
    theta1ALast=theta1ACurr;
    theta1BLast=theta1BCurr;

    equilibrium = !active;
    loopTotal=lpNum;

    std::ofstream outSummary(outDir + "summary.txt");
    const auto tMCEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
    outSummary << "total mc time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    outSummary << "total loop number: " << loopTotal << std::endl;
    outSummary << "nEqCounterStart=" << nEqCounterStart << std::endl;
    outSummary << "equilibrium reached: " << equilibrium << std::endl;
    outSummary << "same: " << same << std::endl;

    outSummary << "lag=" << lag << std::endl;
    outSummary<<"step length="<<stddev<<std::endl;
    outSummary<<"collected number of data points: "<<dataNumInEq<<std::endl;
    outSummary.close();


}


///
/// @param lag
/// @param loopEq
/// @param rInit
/// @param theta0BInit
/// @param theta1AInit
/// @param theta1BInit
/// @param U_ptr
/// @param r_ptr
/// @param theta0B_ptr
/// @param theta1A_ptr
/// @param theta1B_ptr
void version1RingQuadratic::executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const double &rInit,
                                               const double &theta0BInit, const double &theta1AInit, const double &theta1BInit,double * U_ptr,double *r_ptr, double *theta0B_ptr, double *theta1A_ptr, double * theta1B_ptr ) {


    if (dataNumTotal <= dataNumInEq) {
        return;
    }

    unsigned long long remainingDataNum = this->dataNumTotal - this->dataNumInEq;
    unsigned long long remainingLoopNum = remainingDataNum * lag;
    std::cout << "remainingDataNum=" << remainingDataNum << std::endl;

    std::cout << "remainingLoopNum=" << remainingLoopNum << std::endl;
    unsigned long long lastLoopNum = remainingLoopNum % writeInterval;

    std::cout<<"lastLoopNum="<<lastLoopNum<<std::endl;




    double rCurr(rInit);
    double theta0BCurr(theta0BInit);
    double theta1ACurr(theta1AInit);
    double theta1BCurr(theta1BInit);


    double UCurr = (*potFuncPtr)(rCurr,theta0BCurr,theta1ACurr,theta1BCurr);
    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    //output directory
    std::ostringstream sObjT;
    sObjT << std::fixed;
    sObjT << std::setprecision(10);
    sObjT << T;
    std::string TStr = sObjT.str();
    std::string funcName = demangle(typeid(*potFuncPtr).name());

//    std::string  initFuncName= demangle(typeid(initFuncName).name());
    std::string outDir = "./version1Data/1d/func" + funcName + "/row" + std::to_string(rowNum) + "/T" + TStr + "/";
    std::string outUAllPickleSubDir = outDir + "UAllPickle/";
    std::string outUAllBinSubDir = outDir + "UAllBin/";

    std::string out_r_AllBinSubDir = outDir + "r_AllBin/";
    std::string out_theta0B_AllBinSubDir=outDir+"theta0B_AllBin/";

    std::string out_theta1A_AllBinSubDir=outDir+"theta1A_AllBin/";

    std::string out_theta1B_AllBinSubDir=outDir+"theta1B_AllBin/";

    const auto tMCStart{std::chrono::steady_clock::now()};
//remainingLoopNum-lastLoopNum
    unsigned long long lpNum = 0;
    while (lpNum < remainingLoopNum - lastLoopNum) {


        //propose a move
        double rNext;
        double theta0BNext;
        double theta1ANext;
        double theta1BNext;
        proposal(rCurr,theta0BCurr,theta1ACurr,theta1BCurr,rNext,theta0BNext,theta1ANext,theta1BNext);
        double r = acceptanceRatio(rCurr,theta0BCurr,theta1ACurr,theta1BCurr,rNext,theta0BNext,theta1ANext,theta1BNext);
        double u = distUnif01(e2);
        if (u <= r) {
            rCurr=rNext;
            theta0BCurr=theta0BNext;
            theta1ACurr=theta1ANext;
            theta1BCurr=theta1BNext;
            UCurr = (*potFuncPtr)(rCurr,theta0BCurr,theta1ACurr,theta1BCurr);
        }//end of accept-reject

        unsigned long long indOfArrayU = lpNum % writeInterval;


        U_ptr[indOfArrayU] = UCurr;
        r_ptr[indOfArrayU]=rCurr;
        theta0B_ptr[indOfArrayU]=theta0BCurr;
        theta1A_ptr[indOfArrayU]=theta1ACurr;
        theta1B_ptr[indOfArrayU]=theta1BCurr;




        //write to file every writeInterval loops

        if ((lpNum + 1) % loopToWrite == 0 and lpNum > 1) {
            unsigned long long sizeOfArrayU = indOfArrayU + 1;

            std::string filenameMiddle = "loopEnd" + std::to_string(lpNum);
            std::string outUPicleFileName = outUAllPickleSubDir + filenameMiddle + ".UAll.pkl";
            std::string outUBinFileName = outUAllBinSubDir + filenameMiddle + ".UAll.bin";
            save_array_to_pickle(U_ptr, sizeOfArrayU, outUPicleFileName);
            save_to_bin_file(U_ptr, sizeOfArrayU, outUBinFileName);

            std::string out_r_BinFileName=out_r_AllBinSubDir+filenameMiddle+".r.bin";
            save_to_bin_file(r_ptr,sizeOfArrayU,out_r_BinFileName);

            std::string out_theta0B_BinFileName=out_theta0B_AllBinSubDir+filenameMiddle+".theta0B_All.bin";
            save_to_bin_file(theta0B_ptr,sizeOfArrayU,out_theta0B_BinFileName);

            std::string out_theta1A_BinFileName=out_theta1A_AllBinSubDir+filenameMiddle+".theta1A_All.bin";
            save_to_bin_file(theta1A_ptr,sizeOfArrayU,out_theta1A_BinFileName);

            std::string out_theta1B_BinFileName=out_theta1B_AllBinSubDir+filenameMiddle+".theta1B_All.bin";
            save_to_bin_file(theta1B_ptr,sizeOfArrayU,out_theta1B_BinFileName);

            const auto tWriteEnd{std::chrono::steady_clock::now()};

            const std::chrono::duration<double> elapsed_seconds{tWriteEnd - tMCStart};
            std::cout << "loop " << lpNum << std::endl;
            std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;


        }//end write (if)

        lpNum++;
    }//end while

    for(unsigned long long lpFinal=0;lpFinal<lastLoopNum;lpFinal++){

        //propose a move
        double rNext;
        double theta0BNext;
        double theta1ANext;
        double theta1BNext;
        proposal(rCurr,theta0BCurr,theta1ACurr,theta1BCurr,rNext,theta0BNext,theta1ANext,theta1BNext);
        double r = acceptanceRatio(rCurr,theta0BCurr,theta1ACurr,theta1BCurr,rNext,theta0BNext,theta1ANext,theta1BNext);
        double u = distUnif01(e2);
        if (u <= r) {
            rCurr=rNext;
            theta0BCurr=theta0BNext;
            theta1ACurr=theta1ANext;
            theta1BCurr=theta1BNext;
            UCurr = (*potFuncPtr)(rCurr,theta0BCurr,theta1ACurr,theta1BCurr);
        }//end of accept-reject


        U_ptr[lpFinal]=UCurr;
        r_ptr[lpFinal]=rCurr;
        theta0B_ptr[lpFinal]=theta0BCurr;
        theta1A_ptr[lpFinal]=theta1ACurr;
        theta1B_ptr[lpFinal]=theta1BCurr;




    }//end of final for loops

    std::string filenameMiddle = "loopEnd" + std::to_string(lpNum+lastLoopNum-1);
    std::string outUPicleFileName = outUAllPickleSubDir + filenameMiddle + ".UAll.pkl";
    std::string outUBinFileName = outUAllBinSubDir + filenameMiddle + ".UAll.bin";

    save_array_to_pickle(U_ptr,lastLoopNum,outUPicleFileName);
    save_to_bin_file(U_ptr,lastLoopNum,outUBinFileName);

    std::string out_theta0B_BinFileName=out_theta0B_AllBinSubDir+filenameMiddle+".theta0B_All.bin";
    save_to_bin_file(theta0B_ptr,lastLoopNum,out_theta0B_BinFileName);

    std::string out_theta1A_BinFileName=out_theta1A_AllBinSubDir+filenameMiddle+".theta1A_All.bin";
    save_to_bin_file(theta1A_ptr,lastLoopNum,out_theta1A_BinFileName);

    std::string out_theta1B_BinFileName=out_theta1B_AllBinSubDir+filenameMiddle+".theta1B_All.bin";
    save_to_bin_file(theta1B_ptr,lastLoopNum,out_theta1B_BinFileName);

    const auto tMCEnd{std::chrono::steady_clock::now()};

    const std::chrono::duration<double> elapsed_seconds{tMCEnd - tMCStart};
    std::cout << "end " << lpNum << std::endl;
    std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;



}//end of function executionMCAfterEq


void save_array_to_pickle(double *ptr, std::size_t size, const std::string& filename) {
    using namespace boost::python;
    try {
        Py_Initialize();  // Initialize the Python interpreter
        if (!Py_IsInitialized()) {
            throw std::runtime_error("Failed to initialize Python interpreter");
        }

        // Debug output
        std::cout << "Python interpreter initialized successfully." << std::endl;

        // Import the pickle module
        object pickle = import("pickle");
        object pickle_dumps = pickle.attr("dumps");

        // Create a Python list from the C++ array
        list py_list;
        for (std::size_t i = 0; i < size; ++i) {
            py_list.append(ptr[i]);
        }

        // Serialize the list using pickle.dumps
        object serialized_array = pickle_dumps(py_list);

        // Extract the serialized data as a string
        std::string serialized_str = extract<std::string>(serialized_array);

        // Write the serialized data to a file
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open file for writing");
        }
        file.write(serialized_str.data(), serialized_str.size());
        file.close();

        // Debug output
        std::cout << "Array serialized and written to file successfully." << std::endl;
    } catch (const error_already_set&) {
        PyErr_Print();
        std::cerr << "Boost.Python error occurred." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    if (Py_IsInitialized()) {
        Py_Finalize();  // Finalize the Python interpreter
    }
}


///to msgpack bin file
void save_to_bin_file(double * data, unsigned long long  size, const std::string& filename){


/// Create a MessagePack buffer
    msgpack::sbuffer sbuf;

    // Pack the array of doubles
    msgpack::packer<msgpack::sbuffer> packer(sbuf);
    packer.pack_array(size);
    for (unsigned long long i = 0; i < size; ++i) {
        packer.pack_double(data[i]);
    }

    // Write the packed data to a file
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outfile.write(sbuf.data(), sbuf.size());
    outfile.close();

}