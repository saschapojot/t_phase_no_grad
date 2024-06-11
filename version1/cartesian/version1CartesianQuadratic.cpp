//
// Created by polya on 5/23/24.
//
#include "version1CartesianQuadratic.hpp"


void version1CartesianQuadratic::parseCSV(const int &rowNum, double &a1, double &a2, double &c1, double &c2,double &mA, double &mB) {

    std::string filePath = "./version1Input/1d/cartesian/cartesianQuadratic.csv";
    std::string pyFile = "./version1/cartesian/readCSV.py";
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



    std::cout<<"a1="<<a1<<", a2="<<a2<<", c1="<<c1<<", c2="<<c2<<std::endl;






//std::cout<<"r0="<<r0<<std::endl;



}

///
/// @param cmd python execution string
/// @return signal from the python
std::string version1CartesianQuadratic::execPython(const char *cmd) {
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
/// @param L
/// @param x0A
/// @param x0B
/// @param x1A
/// @param x1B
/// @return
double version1CartesianQuadratic::f(const double &L,const double& x0A, const double &x0B, const double&x1A, const double&x1B){
    return this->beta * ((*potFuncPtr)(L,x0A ,x0B,x1A,x1B));



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
void version1CartesianQuadratic::proposal(const double &LCurr, const double& x0ACurr,const double&x0BCurr, const double& x1ACurr,const double &x1BCurr,
                                          double & LNext, double &x0ANext, double &x0BNext, double &x1ANext, double &x1BNext) {


    //next L
    LNext=generate_nearby_positive_value(LCurr,stddev);

    //x

    //x0A

    x0ANext=generate_nearby_0_L(x0ACurr,stddev,LNext);
    //x0B

    x0BNext=generate_nearby_0_L(x0BCurr,stddev,LNext);

    //x1A
    x1ANext= generate_nearby_0_L(x1ACurr,stddev,LNext);

    //x1B
    x1BNext= generate_nearby_0_L(x1BCurr,stddev,LNext);




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
/// @return
double version1CartesianQuadratic::acceptanceRatio(const double &LCurr,const double &x0ACurr, const double &x0BCurr, const double& x1ACurr, const double &x1BCurr,
                                                   const double &LNext, const double&x0ANext, const double &x0BNext, const double &x1ANext, const double &x1BNext){


    double numerator = -f(LNext,x0ANext, x0BNext, x1ANext, x1BNext);

    double denominator = -f(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
//    double UCurr=(*potFuncPtr)(xA, xB,LCurr);
//    double UNext=(*potFuncPtr)(zA, zB, LNext);
//    std::cout<<"UCurr="<<UCurr<<", UNext="<<UNext<<std::endl;
    double ratio = std::exp(numerator - denominator);

    return std::min(1.0, ratio);

}



///
/// @param LInit
/// @param x0AInit
/// @param x0BInit
/// @param x1AInit
/// @param x1BInit
void version1CartesianQuadratic::initPositionsEquiDistance(double &LInit,double &x0AInit, double &x0BInit,  double &x1AInit, double &x1BInit){

    double a = 10;
    LInit=2*a;

    x0AInit=0.1*LInit;

    x0BInit=0.2*LInit;

    x1AInit=0.5*LInit;

    x1BInit=0.9*LInit;



}


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
void version1CartesianQuadratic::readEqMc(unsigned long long &lag,  unsigned long long &loopTotal,bool &equilibrium, bool& same, double &LLast,
                                          double &x0ALast, double &x0BLast,  double &x1ALast, double &x1BLast,double * U_ptr,double *L_ptr, double *x0A_ptr,double *x0B_ptr, double *x1A_ptr, double * x1B_ptr){

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

    std::string out_L_AllBinSubDir = outDir + "L_AllBin/";
    std::string out_L_AllPickleSubDir=outDir+"L_AllPickle/";

    std::string out_x0A_AllBinSubDir=outDir+"x0A_AllBin/";
    std::string out_x0A_AllPickleSubDir=outDir+"x0A_AllPickle/";

    std::string out_x0B_AllBinSubDir=outDir+"x0B_AllBin/";
    std::string out_x0B_AllPickleSubDir=outDir+"x0B_AllPickle/";

    std::string out_x1A_AllBinSubDir=outDir+"x1A_AllBin/";
    std::string out_x1A_AllPickleSubDir=outDir+"x1A_AllPickle/";

    std::string out_x1B_AllBinSubDir=outDir+"x1B_AllBin/";
    std::string out_x1B_AllPickleSubDir=outDir+"x1B_AllPickle/";



    if (!fs::is_directory(outUAllPickleSubDir) || !fs::exists(outUAllPickleSubDir)) {
        fs::create_directories(outUAllPickleSubDir);
    }
    if (!fs::is_directory(outUAllBinSubDir) || !fs::exists(outUAllBinSubDir)) {
        fs::create_directories(outUAllBinSubDir);
    }
    if (!fs::is_directory(out_L_AllBinSubDir) || !fs::exists(out_L_AllBinSubDir)) {
        fs::create_directories(out_L_AllBinSubDir);
    }

    if (!fs::is_directory(out_L_AllPickleSubDir) || !fs::exists(out_L_AllPickleSubDir)) {
        fs::create_directories(out_L_AllPickleSubDir);
    }

    if (!fs::is_directory(out_x0A_AllBinSubDir) || !fs::exists(out_x0A_AllBinSubDir)) {
        fs::create_directories(out_x0A_AllBinSubDir);
    }



    if (!fs::is_directory(out_x0B_AllBinSubDir ) || !fs::exists(out_x0B_AllBinSubDir )) {
        fs::create_directories(out_x0B_AllBinSubDir );
    }

    if (!fs::is_directory(out_x1A_AllBinSubDir ) || !fs::exists(out_x1A_AllBinSubDir )) {
        fs::create_directories(out_x1A_AllBinSubDir );
    }

    if (!fs::is_directory(out_x1B_AllBinSubDir ) || !fs::exists(out_x1B_AllBinSubDir )) {
        fs::create_directories(out_x1B_AllBinSubDir );
    }

    //pickle
    if (!fs::is_directory(out_x0A_AllPickleSubDir ) || !fs::exists(out_x0A_AllPickleSubDir )) {
        fs::create_directories(out_x0A_AllPickleSubDir );
    }
    if (!fs::is_directory(out_x0B_AllPickleSubDir ) || !fs::exists(out_x0B_AllPickleSubDir )) {
        fs::create_directories(out_x0B_AllPickleSubDir );
    }
    if (!fs::is_directory(out_x1A_AllPickleSubDir ) || !fs::exists(out_x1A_AllPickleSubDir )) {
        fs::create_directories(out_x1A_AllPickleSubDir );
    }
    if (!fs::is_directory(out_x1B_AllPickleSubDir ) || !fs::exists(out_x1B_AllPickleSubDir )) {
        fs::create_directories(out_x1B_AllPickleSubDir );
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
    std::smatch match_L_Lag;
    std::smatch match_x0A_Lag;
    std::smatch match_x0B_Lag;
    std::smatch  match_x1A_Lag;
    std::smatch  match_x1B_Lag;

    std::smatch matchUSame;

    std::smatch matchUEq;
    std::smatch match_L_Eq;
    std::smatch match_x0A_Eq;
    std::smatch  match_x0B_Eq;
    std::smatch match_x1A_Eq;
    std::smatch match_x1B_Eq;

    std::smatch matchUCounterStart;
    std::smatch match_L_CounterStart;
    std::smatch match_x0A_CounterStart;
    std::smatch match_x0B_CounterStart;
    std::smatch match_x1A_CounterStart;
    std::smatch match_x1B_CounterStart;

    std::smatch matchUDataNumEq;
    std::smatch match_L_DataNumEq;
    std::smatch match_x0A_DataNumEq;
    std::smatch match_x0B_DataNumEq;
    std::smatch match_x1A_DataNumEq;
    std::smatch match_x1B_DataNumEq;








    double LCurr;
    double x0ACurr;
    double x0BCurr;
    double x1ACurr;
    double x1BCurr;


    this->initPositionsEquiDistance(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
    double UCurr = (*potFuncPtr)(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
    int lpNum=0;
    bool active = true;
    const auto tMCStart{std::chrono::steady_clock::now()};



    while(lpNum<this->loopMax and active==true){
        //propose a move
        double LNext;
        double x0ANext;
        double x0BNext;
        double x1ANext;
        double x1BNext;

        proposal(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,LNext,x0ANext,x0BNext,x1ANext,x1BNext);
        double r = acceptanceRatio(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,LNext,x0ANext,x0BNext,x1ANext,x1BNext);
        double u = distUnif01(e2);
//        double UTmp=UCurr;
        if (u <= r) {
//                std::cout<<"UCurr="<<UCurr<<std::endl;
            LCurr=LNext;
            x0ACurr=x0ANext;
            x0BCurr=x0BNext;
            x1ACurr=x1ANext;
            x1BCurr=x1BNext;
            UCurr = (*potFuncPtr)(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
//                std::cout<<"UNext="<<UCurr<<std::endl;
//            double UDiff=UCurr-UTmp;
//            std::cout<<"UDiff="<<UDiff<<std::endl;
//            std::cout<<"-UDiff/T="<<-UDiff/T<<std::endl;
        }//end of accept-reject

        U_ptr[lpNum]=UCurr;
        L_ptr[lpNum]=LCurr;
        x0A_ptr[lpNum]=x0ACurr;
        x0B_ptr[lpNum]=x0BCurr;
        x1A_ptr[lpNum]=x1ACurr;
        x1B_ptr[lpNum]=x1BCurr;

        //write to file every loopToWrite loops, and inquire equilibrium
        if ((lpNum+1)%loopToWrite==0 and lpNum>1){
            unsigned long long sizeOfArray=lpNum+1;

//            unsigned long long lpStart=0;

            std::string filenameMiddle="loopStart0ReachEq";
            //U
            std::string outUPickleFileName=outUAllPickleSubDir+filenameMiddle+ ".UAll.pkl";
            std::string outUBinFileName=outUAllBinSubDir+filenameMiddle+".UAll.bin";
            save_array_to_pickle(U_ptr,sizeOfArray,outUPickleFileName);
            save_to_bin_file(U_ptr,sizeOfArray,outUBinFileName);

            //L
            std::string out_L_BinFileName=out_L_AllBinSubDir+filenameMiddle+".L.bin";
            save_to_bin_file(L_ptr,sizeOfArray,out_L_BinFileName);

            std::string out_L_PickleFileName=out_L_AllPickleSubDir+filenameMiddle+".L.pkl";
            save_array_to_pickle(L_ptr,sizeOfArray,out_L_PickleFileName);

            //x0A
            std::string out_x0A_BinFileName=out_x0A_AllBinSubDir+filenameMiddle+".x0A_All.bin";
            save_to_bin_file(x0A_ptr,sizeOfArray,out_x0A_BinFileName);

            std::string out_x0A_PickleFileName=out_x0A_AllPickleSubDir+filenameMiddle+".x0A_All.pkl";
            save_array_to_pickle(x0A_ptr,sizeOfArray,out_x0A_PickleFileName);

            //x0B
            std::string out_x0B_BinFileName=out_x0B_AllBinSubDir+filenameMiddle+".x0B_All.bin";
            save_to_bin_file(x0B_ptr,sizeOfArray,out_x0B_BinFileName);

            std::string out_x0B_PickleFileName=out_x0B_AllPickleSubDir+filenameMiddle+".x0B_All.pkl";
            save_array_to_pickle(x0B_ptr,sizeOfArray,out_x0B_PickleFileName);

            //x1A
            std::string out_x1A_BinFileName=out_x1A_AllBinSubDir+filenameMiddle+".x1A_All.bin";
            save_to_bin_file(x1A_ptr,sizeOfArray,out_x1A_BinFileName);

            std::string out_x1A_PickleFileName=out_x1A_AllPickleSubDir+filenameMiddle+".x1A_All.pkl";
            save_array_to_pickle(x1A_ptr,sizeOfArray,out_x1A_PickleFileName);

            //x1B
            std::string out_x1B_BinFileName=out_x1B_AllBinSubDir+filenameMiddle+".x1B_All.bin";
            save_to_bin_file(x1B_ptr,sizeOfArray,out_x1B_BinFileName);

            std::string out_x1B_PickleFileName=out_x1B_AllPickleSubDir+filenameMiddle+".x1B_All.pkl";
            save_array_to_pickle(x1B_ptr,sizeOfArray,out_x1B_PickleFileName);

            const auto tWriteEnd{std::chrono::steady_clock::now()};

            const std::chrono::duration<double> elapsed_seconds{tWriteEnd - tMCStart};
            std::cout << "loop " << lpNum << std::endl;
            std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;

            //communicate with python to inquire equilibrium

            //inquire equilibrium of U
            std::string commandU = "python3 checkVec.py " + outUPickleFileName;
            std::string resultU;

            //inquire L
            std::string commandL="python3 checkVec.py "+out_L_PickleFileName;
            std::string resultL;

            //inquire x0A
            std::string command_x0A="python3 checkVec.py "+ out_x0A_PickleFileName;
            std::string result_x0A;

            //inquire x0B
            std::string command_x0B="python3 checkVec.py "+out_x0B_PickleFileName;
            std::string result_x0B;

            //inquire x1A
            std::string command_x1A="python3 checkVec.py "+out_x1A_PickleFileName;
            std::string result_x1A;

            //inquire x1B
            std::string command_x1B="python3 checkVec.py "+out_x1B_PickleFileName;
            std::string result_x1B;

            try {
                const auto tPyStart{std::chrono::steady_clock::now()};
                resultU = this->execPython(commandU.c_str());
                std::cout << "U message from python: " << resultU << std::endl;

                resultL=this->execPython(commandL.c_str());
                std::cout<<"L message from python: "<<resultL<<std::endl;

                result_x0A=this->execPython(command_x0A.c_str());
                std::cout<<"x0A message from python: "<<result_x0A<<std::endl;

                result_x0B=this->execPython(command_x0B.c_str());
                std::cout<<"x0B message from python: "<<result_x0B<<std::endl;

                result_x1A=this->execPython(command_x1A.c_str());
                std::cout<<"x1A message from python: "<<result_x1A<<std::endl;

                result_x1B=this->execPython(command_x1B.c_str());
                std::cout<<"x1B message from python: "<<result_x1B<<std::endl;

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


            if (
                    std::regex_search(resultU, matchUEq, eqRegex)
            and std::regex_search(resultL,match_L_Eq,eqRegex)
            and std::regex_search(result_x0A,match_x0A_Eq,eqRegex)
            and std::regex_search(result_x0B,match_x0B_Eq,eqRegex)
            and std::regex_search(result_x1A,match_x1A_Eq,eqRegex)
            and std::regex_search(result_x1B,match_x1B_Eq,eqRegex)

            ) {
                unsigned long long lagU;
                unsigned long long lag_L;
                unsigned long long lag_x0A;
                unsigned long long lag_x0B;
                unsigned long long lag_x1A;
                unsigned long long lag_x1B;

                unsigned long long UnEqCounterStart;
                unsigned long long L_nEqCounterStart;
                unsigned long long x0A_nEqCounterStart;
                unsigned long long x0B_nEqCounterStart;
                unsigned long long x1A_nEqCounterStart;
                unsigned long long x1B_nEqCounterStart;


                unsigned long long  UdataNumInEq;
                 unsigned long long  L_dataNumInEq;
                unsigned long long  x0A_dataNumInEq;
                unsigned long long  x0B_dataNumInEq;
                unsigned long long  x1A_dataNumInEq;
                unsigned long long  x1B_dataNumInEq;

                /////////////////////////////////////////////////////
                //parse U
                //match U lag
                if (std::regex_search(resultU, matchULag, lagRegex)) {


                    std::string lagStrU = matchULag.str(1);
                     lagU = std::stoull(lagStrU);
                    std::cout << "lagU=" << lagU << std::endl;


//

                }

                //match UCounterStart
                if( std::regex_search(resultU,matchUCounterStart,ctStartRegex)){
                    UnEqCounterStart=std::stoull(matchUCounterStart.str(1));

                    std::cout<<"UnEqCounterStart="<<UnEqCounterStart<<std::endl;

                }

                //match UdataNumInEq
                if(std::regex_search(resultU,matchUDataNumEq,dataNumEqRegex)){
                    UdataNumInEq=std::stoull(matchUDataNumEq.str(1));

                    std::cout<<"UdataNumInEq="<<UdataNumInEq<<std::endl;

                }
                // end of parse U
                ///////////////////////////////////////////////////////
                //////////////////////////////////////////////////////
                //parse L

                //match L lag
                if(std::regex_search(resultL,match_L_Lag,lagRegex)){
                std::string lagStrL=match_L_Lag.str(1);
                    lag_L=std::stoull(lagStrL);
                    std::cout<<"lag_L="<<lag_L<<std::endl;


                }

                //match L_nEqCounterStart
                if(std::regex_search(resultL,match_L_CounterStart,ctStartRegex)){
                    L_nEqCounterStart=std::stoull(match_L_CounterStart.str(1));
                    std::cout<<"L_nEqCounterStart="<<L_nEqCounterStart<<std::endl;
                }

                //match L_dataNumInEq
                if(std::regex_search(resultL,match_L_DataNumEq,dataNumEqRegex)){
                    L_dataNumInEq=std::stoull(match_L_DataNumEq.str(1));
                    std::cout<<"L_dataNumInEq="<<L_dataNumInEq<<std::endl;
                }

                // end of parse L
                //////////////////////////////////////

                /////////////////////////////////////
                //parse x0A

                //match x0A lag
                if(std::regex_search(result_x0A,match_x0A_Lag,lagRegex)){
                    std::string lagStr_x0A=match_x0A_Lag.str(1);
                    lag_x0A=std::stoull(lagStr_x0A);
                    std::cout<<"lag_x0A="<<lag_x0A<<std::endl;


                }

                //match x0A_nEqCounterStart
                if(std::regex_search(result_x0A,match_x0A_CounterStart,ctStartRegex)){
                    std::string x0A_ctStartStr=match_x0A_CounterStart.str(1);
                    x0A_nEqCounterStart=std::stoull(x0A_ctStartStr);
                    std::cout<<"x0A_nEqCounterStart="<<x0A_nEqCounterStart<<std::endl;

                }

                // match x0A_dataNumInEq

                if(std::regex_search(result_x0A,match_x0A_DataNumEq,dataNumEqRegex)){
                    x0A_dataNumInEq=std::stoull(match_x0A_DataNumEq.str(1));
                    std::cout<<"x0A_dataNumInEq="<<x0A_dataNumInEq<<std::endl;
                }

                // end of parse x0A
                ///////////////////////////////////

                ///////////////////////////////////
                //parse x0B

                //match x0B lag
                if(std::regex_search(result_x0B,match_x0B_Lag,lagRegex)){
                    std::string lagStr_x0B=match_x0B_Lag.str(1);
                    lag_x0B=std::stoull(lagStr_x0B);
                    std::cout<<"lag_x0B="<<lag_x0B<<std::endl;

                }

                //match x0B_nEqCounterStart
                if(std::regex_search(result_x0B,match_x0B_CounterStart,ctStartRegex)){
                    std::string x0B_ctStartStr=match_x0B_CounterStart.str(1);
                    x0B_nEqCounterStart=std::stoull(x0B_ctStartStr);
                    std::cout<<"x0B_nEqCounterStart="<<x0B_nEqCounterStart<<std::endl;

                }

                // match x0B_dataNumInEq
                if(std::regex_search(result_x0B,match_x0B_DataNumEq,dataNumEqRegex)){

                    x0B_dataNumInEq=std::stoull(match_x0B_DataNumEq.str(1));
                    std::cout<<"x0B_dataNumInEq="<<x0B_dataNumInEq<<std::endl;
                }
                // end of parse x0B
                //////////////////////////////////

                //////////////////////////////////
                //parse x1A

                //match x1A lag
                if(std::regex_search(result_x1A,match_x1A_Lag,lagRegex)){
                    std::string lagStr_x1A=match_x1A_Lag.str(1);
                    lag_x1A=std::stoull(lagStr_x1A);
                    std::cout<<"lag_x1A="<<lag_x1A<<std::endl;

                }

                //match x1A_nEqCounterStart
                if(std::regex_search(result_x1A,match_x1A_CounterStart,ctStartRegex)){
                    std::string x1A_ctStartStr=match_x1A_CounterStart.str(1);
                    x1A_nEqCounterStart=std::stoull(x1A_ctStartStr);
                    std::cout<<"x1A_nEqCounterStart="<<x1A_nEqCounterStart<<std::endl;


                }

                // match x1A_dataNumInEq
                if(std::regex_search(result_x1A,match_x1A_DataNumEq,dataNumEqRegex)){
                    x1A_dataNumInEq=std::stoull(match_x1A_DataNumEq.str(1));
                    std::cout<<"x1A_dataNumInEq="<<x1A_dataNumInEq<<std::endl;

                }
                // end of parse x1A
                //////////////////////////////////

                ///////////////////////////////////
                //parse x1B

                //match x1B lag
                if(std::regex_search(result_x1B,match_x1B_Lag,lagRegex)){
                    std::string lagStr_x1B=match_x1B_Lag.str(1);
                    lag_x1B=std::stoull(lagStr_x1B);
                    std::cout<<"lag_x1B="<<lag_x1B<<std::endl;

                }

                //match x1B_nEqCounterStart
                if(std::regex_search(result_x1B,match_x1B_CounterStart,ctStartRegex)){
                    std::string x1B_ctStartStr=match_x1B_CounterStart.str(1);
                    x1B_nEqCounterStart=std::stoull(x1B_ctStartStr);


                }

                // match x1B_dataNumInEq
                if(std::regex_search(result_x1B,match_x1B_DataNumEq,dataNumEqRegex)){
                    x1B_dataNumInEq=std::stoull(match_x1B_DataNumEq.str(1));
                    std::cout<<"x1B_dataNumInEq="<<x1B_dataNumInEq<<std::endl;
                }

                // end of parse x1B
                //////////////////////////////////


                //use the max of lagU, lag_L, lag_x0A, lag_x0B, lag_x1A, lag_x1B

                std::initializer_list<unsigned long  long > lagsAll=
                        {lagU,lag_L,lag_x0A,lag_x0B,lag_x1A,lag_x1B};


                unsigned long  long maxLag=*std::max_element(lagsAll.begin(),lagsAll.end());

                lag=maxLag;
                std::cout<<"lag is chosen as "<<lag<<std::endl;

                std::initializer_list<unsigned long long > nCounterStartAll=
                        {UnEqCounterStart,L_nEqCounterStart,x0A_nEqCounterStart,x0B_nEqCounterStart,x1A_nEqCounterStart,x1B_nEqCounterStart};
                unsigned long long nCounterStartMax=*std::max_element(nCounterStartAll.begin(),nCounterStartAll.end());

                nEqCounterStart=nCounterStartMax;
                std::cout<<"nEqCounterStart is chosen as "<<nEqCounterStart<<std::endl;

                std::initializer_list<unsigned long long > numDataPointsAll=
                        {UdataNumInEq,L_dataNumInEq,x0A_dataNumInEq,x0B_dataNumInEq,x1A_dataNumInEq,x1B_dataNumInEq};
                unsigned long long  numDataPointsMin=*std::min_element(numDataPointsAll.begin(),numDataPointsAll.end());

                dataNumInEq=numDataPointsMin;

                std::cout<<"dataNumInEq is chosen as "<<dataNumInEq<<std::endl;


                active = false;

            }//end of regex search for equilibrium






        }//end write and inquire (if)




        lpNum++;

    }//end while

    LLast=LCurr;
    x0ALast=x0ACurr;
    x0BLast=x0BCurr;
    x1ALast=x1ACurr;
    x1BLast=x1BCurr;

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
void version1CartesianQuadratic::executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const double &LInit,
                                                    const double &x0AInit, const double &x0BInit, const double &x1AInit, const double &x1BInit,double * U_ptr,double *L_ptr,double *x0A_ptr, double *x0B_ptr, double *x1A_ptr, double * x1B_ptr) {


    if (dataNumTotal <= dataNumInEq) {
        return;
    }

    unsigned long long remainingDataNum = this->dataNumTotal - this->dataNumInEq;
    unsigned long long remainingLoopNum = remainingDataNum * lag;
    std::cout << "remainingDataNum=" << remainingDataNum << std::endl;

    std::cout << "remainingLoopNum=" << remainingLoopNum << std::endl;
    unsigned long long lastLoopNum = remainingLoopNum % writeInterval;

    std::cout<<"lastLoopNum="<<lastLoopNum<<std::endl;




    double LCurr(LInit);
    double x0ACurr(x0AInit);
    double x0BCurr(x0BInit);
    double x1ACurr(x1AInit);
    double x1BCurr(x1BInit);


    double UCurr = (*potFuncPtr)(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
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

    std::string out_L_AllBinSubDir = outDir + "L_AllBin/";

    std::string out_x0A_AllBinSubDir=outDir+"x0A_AllBin/";
    std::string out_x0B_AllBinSubDir=outDir+"x0B_AllBin/";

    std::string out_x1A_AllBinSubDir=outDir+"x1A_AllBin/";

    std::string out_x1B_AllBinSubDir=outDir+"x1B_AllBin/";

    const auto tMCStart{std::chrono::steady_clock::now()};
//remainingLoopNum-lastLoopNum
    unsigned long long lpNum = 0;
    while (lpNum < remainingLoopNum - lastLoopNum) {


        //propose a move
        double LNext;
        double x0ANext;
        double x0BNext;
        double x1ANext;
        double x1BNext;
        proposal(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,LNext,x0ANext,x0BNext,x1ANext,x1BNext);
        double r = acceptanceRatio(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,LNext,x0ANext,x0BNext,x1ANext,x1BNext);
        double u = distUnif01(e2);
        if (u <= r) {
            LCurr=LNext;
            x0ACurr=x0ANext;
            x0BCurr=x0BNext;
            x1ACurr=x1ANext;
            x1BCurr=x1BNext;
            UCurr = (*potFuncPtr)(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
        }//end of accept-reject

        unsigned long long indOfArrayU = lpNum % writeInterval;


        U_ptr[indOfArrayU] = UCurr;
        L_ptr[indOfArrayU]=LCurr;
        x0A_ptr[indOfArrayU]=x0ACurr;
        x0B_ptr[indOfArrayU]=x0BCurr;
        x1A_ptr[indOfArrayU]=x1ACurr;
        x1B_ptr[indOfArrayU]=x1BCurr;




        //write to file every writeInterval loops

        if ((lpNum + 1) % loopToWrite == 0 and lpNum > 1) {
            unsigned long long sizeOfArrayU = indOfArrayU + 1;

            std::string filenameMiddle = "loopEnd" + std::to_string(lpNum);
            std::string outUPicleFileName = outUAllPickleSubDir + filenameMiddle + ".UAll.pkl";
            std::string outUBinFileName = outUAllBinSubDir + filenameMiddle + ".UAll.bin";
            save_array_to_pickle(U_ptr, sizeOfArrayU, outUPicleFileName);
            save_to_bin_file(U_ptr, sizeOfArrayU, outUBinFileName);

            std::string out_L_BinFileName=out_L_AllBinSubDir+filenameMiddle+".L.bin";
            save_to_bin_file(L_ptr,sizeOfArrayU,out_L_BinFileName);

            std::string out_x0A_BinFileName=out_x0A_AllBinSubDir+filenameMiddle+".x0A_All.bin";
            save_to_bin_file(x0A_ptr,sizeOfArrayU,out_x0A_BinFileName);

            std::string out_x0B_BinFileName=out_x0B_AllBinSubDir+filenameMiddle+".x0B_All.bin";
            save_to_bin_file(x0B_ptr,sizeOfArrayU,out_x0B_BinFileName);

            std::string out_x1A_BinFileName=out_x1A_AllBinSubDir+filenameMiddle+".x1A_All.bin";
            save_to_bin_file(x1A_ptr,sizeOfArrayU,out_x1A_BinFileName);

            std::string out_x1B_BinFileName=out_x1B_AllBinSubDir+filenameMiddle+".x1B_All.bin";
            save_to_bin_file(x1B_ptr,sizeOfArrayU,out_x1B_BinFileName);

            const auto tWriteEnd{std::chrono::steady_clock::now()};

            const std::chrono::duration<double> elapsed_seconds{tWriteEnd - tMCStart};
            std::cout << "loop " << lpNum << std::endl;
            std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;


        }//end write (if)

        lpNum++;
    }//end while

    for(unsigned long long lpFinal=0;lpFinal<lastLoopNum;lpFinal++){

        //propose a move
        double LNext;
        double x0ANext;
        double x0BNext;
        double x1ANext;
        double x1BNext;
        proposal(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,LNext,x0ANext,x0BNext,x1ANext,x1BNext);
        double r = acceptanceRatio(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr,LNext,x0ANext,x0BNext,x1ANext,x1BNext);
        double u = distUnif01(e2);
        if (u <= r) {
            LCurr=LNext;
            x0ACurr=x0ANext;
            x0BCurr=x0BNext;
            x1ACurr=x1ANext;
            x1BCurr=x1BNext;
            UCurr = (*potFuncPtr)(LCurr,x0ACurr,x0BCurr,x1ACurr,x1BCurr);
        }//end of accept-reject


        U_ptr[lpFinal]=UCurr;
        L_ptr[lpFinal]=LCurr;
        x0A_ptr[lpFinal]=x0ACurr;
        x0B_ptr[lpFinal]=x0BCurr;
        x1A_ptr[lpFinal]=x1ACurr;
        x1B_ptr[lpFinal]=x1BCurr;




    }//end of final for loops

    std::string filenameMiddle = "loopEnd" + std::to_string(lpNum+lastLoopNum-1);
    std::string outUPicleFileName = outUAllPickleSubDir + filenameMiddle + ".UAll.pkl";
    std::string outUBinFileName = outUAllBinSubDir + filenameMiddle + ".UAll.bin";

    save_array_to_pickle(U_ptr,lastLoopNum,outUPicleFileName);
    save_to_bin_file(U_ptr,lastLoopNum,outUBinFileName);

    std::string out_L_BinFileName=out_L_AllBinSubDir+filenameMiddle+".L.bin";
    save_to_bin_file(L_ptr,lastLoopNum,out_L_BinFileName);

    std::string out_x0A_BinFileName=out_x0A_AllBinSubDir+filenameMiddle+".x0A_All.bin";
    save_to_bin_file(x0A_ptr,lastLoopNum,out_x0A_BinFileName);

    std::string out_x0B_BinFileName=out_x0B_AllBinSubDir+filenameMiddle+".x0B_All.bin";
    save_to_bin_file(x0B_ptr,lastLoopNum,out_x0B_BinFileName);

    std::string out_x1A_BinFileName=out_x1A_AllBinSubDir+filenameMiddle+".x1A_All.bin";
    save_to_bin_file(x1A_ptr,lastLoopNum,out_x1A_BinFileName);

    std::string out_x1B_BinFileName=out_x1B_AllBinSubDir+filenameMiddle+".x1B_All.bin";
    save_to_bin_file(x1B_ptr,lastLoopNum,out_x1B_BinFileName);

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