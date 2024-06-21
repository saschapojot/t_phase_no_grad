//
// Created by polya on 5/23/24.
//
#include "version1DistRelative.hpp"


void version1DistRelative::parseCSV(const int &rowNum, double &a1, double &a2, double &c1, double &c2,double &mA, double &mB) {

    std::string filePath = "./version1Input/1d/distRelative/distRelative.csv";
    std::string pyFile = "./version1/distRelative/readCSV.py";
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
std::string version1DistRelative::execPython(const char *cmd) {
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
/// @param y0
/// @param z0
/// @param y1
/// @return
double version1DistRelative::f(const double &L,const double& y0, const double &z0, const double&y1){
    return this->beta * ((*potFuncPtr)(L,y0,z0,y1));



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
void version1DistRelative::proposal(const double &LCurr, const double& y0Curr,const double& z0Curr, const double& y1Curr,
                                    double & LNext, double & y0Next, double & z0Next, double & y1Next) {


    //next L
    LNext=generate_nearby_positive_value(LCurr,stddev);
//    LNext = generate_LNext(LCurr, stddev);

   y0Next= generate_nearby_mL_L(y0Curr,stddev,LNext);
//    y0Next = generate_distNext(y0Curr, stddev, LNext);

   z0Next= generate_nearby_mL_L(z0Curr,stddev,LNext);
//    z0Next = generate_distNext(z0Curr, stddev, LNext);

   y1Next= generate_nearby_mL_L(y1Curr,stddev,LNext);
//    y1Next = generate_distNext(y1Curr, stddev, LNext);


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
/// @return
double version1DistRelative::acceptanceRatio(const double &LCurr,const double &y0Curr, const double &z0Curr, const double& y1Curr,
                                             const double &LNext, const double& y0Next, const double & z0Next, const double & y1Next){


    double numerator = -f(LNext, y0Next,z0Next,y1Next);

    double denominator = -f(LCurr, y0Curr,z0Curr,y1Curr);
//    double UCurr=(*potFuncPtr)(xA, xB,LCurr);
//    double UNext=(*potFuncPtr)(zA, zB, LNext);
//    std::cout<<"UCurr="<<UCurr<<", UNext="<<UNext<<std::endl;
    double ratio = std::exp(numerator - denominator);

    return std::min(1.0, ratio);

}



///
/// @param LInit
/// @param y0Init
/// @param z0Init
/// @param y1Init
void version1DistRelative::initPositionsEquiDistance(double &LInit,double &y0Init, double &z0Init,  double &y1Init){

    double a = 10;
    LInit=2*a;

    y0Init=0.1*LInit;

    z0Init=0.2*LInit;

    y1Init=0.5*LInit;





}


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
void version1DistRelative::readEqMc(unsigned long long &lag,  unsigned long long &loopTotal,bool &equilibrium, bool& same, double &LLast,
                                    double &y0Last, double & z0Last,  double & y1Last,double * U_ptr,double *L_ptr, double *y0_ptr,double *z0_ptr, double *y1_ptr){

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

    std::string out_y0_AllBinSubDir=outDir+"y0_AllBin/";
    std::string out_y0_AllPickleSubDir=outDir+"y0_AllPickle/";

    std::string out_z0_AllBinSubDir=outDir+"z0_AllBin/";
    std::string out_z0_AllPickleSubDir=outDir+"z0_AllPickle/";

    std::string out_y1_AllBinSubDir=outDir+"y1_AllBin/";
    std::string out_y1_AllPickleSubDir=outDir+"y1_AllPickle/";





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

    if (!fs::is_directory(out_y0_AllBinSubDir) || !fs::exists(out_y0_AllBinSubDir)) {
        fs::create_directories(out_y0_AllBinSubDir);
    }



    if (!fs::is_directory(out_y0_AllPickleSubDir ) || !fs::exists(out_y0_AllPickleSubDir )) {
        fs::create_directories(out_y0_AllPickleSubDir );
    }

    if (!fs::is_directory(out_z0_AllBinSubDir ) || !fs::exists(out_z0_AllBinSubDir )) {
        fs::create_directories(out_z0_AllBinSubDir );
    }

    if (!fs::is_directory(out_z0_AllPickleSubDir ) || !fs::exists(out_z0_AllPickleSubDir )) {
        fs::create_directories(out_z0_AllPickleSubDir );
    }


    if (!fs::is_directory(out_y1_AllBinSubDir ) || !fs::exists(out_y1_AllBinSubDir )) {
        fs::create_directories(out_y1_AllBinSubDir );
    }
    if (!fs::is_directory(out_y1_AllPickleSubDir ) || !fs::exists(out_y1_AllPickleSubDir )) {
        fs::create_directories(out_y1_AllPickleSubDir );
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
    std::smatch match_y0_Lag;
    std::smatch match_z0_Lag;
    std::smatch  match_y1_Lag;


    std::smatch matchUSame;

    std::smatch matchUEq;
    std::smatch match_L_Eq;
    std::smatch match_y0_Eq;
    std::smatch  match_z0_Eq;
    std::smatch match_y1_Eq;


    std::smatch matchUCounterStart;
    std::smatch match_L_CounterStart;
    std::smatch match_y0_CounterStart;
    std::smatch match_z0_CounterStart;
    std::smatch match_y1_CounterStart;


    std::smatch matchUDataNumEq;
    std::smatch match_L_DataNumEq;
    std::smatch match_y0_DataNumEq;
    std::smatch match_z0_DataNumEq;
    std::smatch match_y1_DataNumEq;









    double LCurr;
    double y0Curr;
    double z0Curr;
    double y1Curr;



    this->initPositionsEquiDistance(LCurr,y0Curr,z0Curr,y1Curr);
    double UCurr = (*potFuncPtr)(LCurr,y0Curr,z0Curr,y1Curr);
    int lpNum=0;
    bool active = true;
    const auto tMCStart{std::chrono::steady_clock::now()};



    while(lpNum<this->loopMax and active==true){
        //propose a move
        double LNext;
        double y0Next;
        double z0Next;
        double y1Next;


        proposal(LCurr,y0Curr, z0Curr,y1Curr,LNext,y0Next,z0Next,y1Next);
        double r = acceptanceRatio(LCurr,y0Curr, z0Curr, y1Curr,LNext, y0Next,z0Next,y1Next);
        double u = distUnif01(e2);
//        double UTmp=UCurr;
        if (u <= r) {
//                std::cout<<"UCurr="<<UCurr<<std::endl;
            LCurr=LNext;
            y0Curr=y0Next;
            z0Curr=z0Next;
            y1Curr=y1Next;

            UCurr = (*potFuncPtr)(LCurr,y0Curr,z0Curr,y1Curr);
//                std::cout<<"UNext="<<UCurr<<std::endl;
//            double UDiff=UCurr-UTmp;
//            std::cout<<"UDiff="<<UDiff<<std::endl;
//            std::cout<<"-UDiff/T="<<-UDiff/T<<std::endl;
        }//end of accept-reject

        U_ptr[lpNum]=UCurr;
        L_ptr[lpNum]=LCurr;
        y0_ptr[lpNum]=y0Curr;
        z0_ptr[lpNum]=z0Curr;
        y1_ptr[lpNum]=y1Curr;


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

            //y0
            std::string out_y0_BinFileName=out_y0_AllBinSubDir+filenameMiddle+".y0_All.bin";
            save_to_bin_file(y0_ptr,sizeOfArray,out_y0_BinFileName);

            std::string out_y0_PickleFileName=out_y0_AllPickleSubDir+filenameMiddle+".y0_All.pkl";
            save_array_to_pickle(y0_ptr,sizeOfArray,out_y0_PickleFileName);

            //z0
            std::string out_z0_BinFileName=out_z0_AllBinSubDir+filenameMiddle+".z0_All.bin";
            save_to_bin_file(z0_ptr,sizeOfArray,out_z0_BinFileName);

            std::string out_z0_PickleFileName=out_z0_AllPickleSubDir+filenameMiddle+".z0_All.pkl";
            save_array_to_pickle(z0_ptr,sizeOfArray,out_z0_PickleFileName);

            //y1
            std::string out_y1_BinFileName=out_y1_AllBinSubDir+filenameMiddle+".y1_All.bin";
            save_to_bin_file(y1_ptr,sizeOfArray,out_y1_BinFileName);

            std::string out_y1_PickleFileName=out_y1_AllPickleSubDir+filenameMiddle+".y1_All.pkl";
            save_array_to_pickle(y1_ptr,sizeOfArray,out_y1_PickleFileName);





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

            //inquire y0
            std::string command_y0="python3 checkVec.py "+ out_y0_PickleFileName;
            std::string result_y0;

            //inquire z0
            std::string command_z0="python3 checkVec.py "+out_z0_PickleFileName;
            std::string result_z0;

            //inquire y1
            std::string command_y1="python3 checkVec.py "+out_y1_PickleFileName;
            std::string result_y1;



            try {
                const auto tPyStart{std::chrono::steady_clock::now()};
                resultU = this->execPython(commandU.c_str());
                std::cout << "U message from python: " << resultU << std::endl;

                resultL=this->execPython(commandL.c_str());
                std::cout<<"L message from python: "<<resultL<<std::endl;

                result_y0=this->execPython(command_y0.c_str());
                std::cout<<"y0 message from python: "<<result_y0<<std::endl;

                result_z0=this->execPython(command_z0.c_str());
                std::cout<<"z0 message from python: "<<result_z0<<std::endl;

                result_y1=this->execPython(command_y1.c_str());
                std::cout<<"y1 message from python: "<<result_y1<<std::endl;



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
            and std::regex_search(result_y0,match_y0_Eq,eqRegex)
            and std::regex_search(result_z0,match_z0_Eq,eqRegex)
            and std::regex_search(result_y1,match_y1_Eq,eqRegex)


            ) {
                unsigned long long lagU;
                unsigned long long lag_L;
                unsigned long long lag_y0;
                unsigned long long lag_z0;
                unsigned long long lag_y1;


                unsigned long long UnEqCounterStart;
                unsigned long long L_nEqCounterStart;
                unsigned long long y0_nEqCounterStart;
                unsigned long long z0_nEqCounterStart;
                unsigned long long y1_nEqCounterStart;



                unsigned long long  UdataNumInEq;
                 unsigned long long  L_dataNumInEq;
                unsigned long long  y0_dataNumInEq;
                unsigned long long  z0_dataNumInEq;
                unsigned long long  y1_dataNumInEq;


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
                //parse y0

                //match y0 lag
                if(std::regex_search(result_y0,match_y0_Lag,lagRegex)){
                    std::string lagStr_y0=match_y0_Lag.str(1);
                    lag_y0=std::stoull(lagStr_y0);
                    std::cout<<"lag_y0="<<lag_y0<<std::endl;


                }

                //match y0_nEqCounterStart
                if(std::regex_search(result_y0,match_y0_CounterStart,ctStartRegex)){
                    std::string y0_ctStartStr=match_y0_CounterStart.str(1);
                    y0_nEqCounterStart=std::stoull(y0_ctStartStr);
                    std::cout<<"y0_nEqCounterStart="<<y0_nEqCounterStart<<std::endl;

                }

                // match y0_dataNumInEq

                if(std::regex_search(result_y0,match_y0_DataNumEq,dataNumEqRegex)){
                    y0_dataNumInEq=std::stoull(match_y0_DataNumEq.str(1));
                    std::cout<<"y0_dataNumInEq="<<y0_dataNumInEq<<std::endl;
                }

                // end of parse y0
                ///////////////////////////////////

                ///////////////////////////////////
                //parse z0

                //match z0 lag
                if(std::regex_search(result_z0,match_z0_Lag,lagRegex)){
                    std::string lagStr_z0=match_z0_Lag.str(1);
                    lag_z0=std::stoull(lagStr_z0);
                    std::cout<<"lag_z0="<<lag_z0<<std::endl;

                }

                //match z0_nEqCounterStart
                if(std::regex_search(result_z0,match_z0_CounterStart,ctStartRegex)){
                    std::string z0_ctStartStr=match_z0_CounterStart.str(1);
                    z0_nEqCounterStart=std::stoull(z0_ctStartStr);
                    std::cout<<"z0_nEqCounterStart="<<z0_nEqCounterStart<<std::endl;

                }

                // match z0_dataNumInEq
                if(std::regex_search(result_z0,match_z0_DataNumEq,dataNumEqRegex)){

                    z0_dataNumInEq=std::stoull(match_z0_DataNumEq.str(1));
                    std::cout<<"z0_dataNumInEq="<<z0_dataNumInEq<<std::endl;
                }
                // end of parse z0
                //////////////////////////////////

                //////////////////////////////////
                //parse y1

                //match  y1 lag
                if(std::regex_search(result_y1,match_y1_Lag,lagRegex)){
                    std::string lagStr_y1=match_y1_Lag.str(1);
                    lag_y1=std::stoull(lagStr_y1);
                    std::cout<<"lag_y1="<<lag_y1<<std::endl;

                }

                //match y1_nEqCounterStart
                if(std::regex_search(result_y1,match_y1_CounterStart,ctStartRegex)){
                    std::string y1_ctStartStr=match_y1_CounterStart.str(1);
                    y1_nEqCounterStart=std::stoull(y1_ctStartStr);
                    std::cout<<"y1_nEqCounterStart="<<y1_nEqCounterStart<<std::endl;


                }

                // match y1_dataNumInEq
                if(std::regex_search(result_y1,match_y1_DataNumEq,dataNumEqRegex)){
                    y1_dataNumInEq=std::stoull(match_y1_DataNumEq.str(1));
                    std::cout<<"y1_dataNumInEq="<<y1_dataNumInEq<<std::endl;

                }
                // end of parse y1
                //////////////////////////////////




                //use the max of lagU, lag_L, lag_y0, lag_z0, lag_y1

                std::initializer_list<unsigned long  long > lagsAll=
                        {lagU,lag_L,lag_y0,lag_z0,lag_y1};


                unsigned long  long maxLag=*std::max_element(lagsAll.begin(),lagsAll.end());

                lag=maxLag;
                std::cout<<"lag is chosen as "<<lag<<std::endl;

                std::initializer_list<unsigned long long > nCounterStartAll=
                        {UnEqCounterStart,L_nEqCounterStart,y0_nEqCounterStart,z0_nEqCounterStart,y1_nEqCounterStart};
                unsigned long long nCounterStartMax=*std::max_element(nCounterStartAll.begin(),nCounterStartAll.end());

                nEqCounterStart=nCounterStartMax;
                std::cout<<"nEqCounterStart is chosen as "<<nEqCounterStart<<std::endl;

                std::initializer_list<unsigned long long > numDataPointsAll=
                        {UdataNumInEq,L_dataNumInEq,y0_dataNumInEq,z0_dataNumInEq,y1_dataNumInEq};
                unsigned long long  numDataPointsMin=*std::min_element(numDataPointsAll.begin(),numDataPointsAll.end());

                dataNumInEq=numDataPointsMin;

                std::cout<<"dataNumInEq is chosen as "<<dataNumInEq<<std::endl;


                active = false;

            }//end of regex search for equilibrium






        }//end write and inquire (if)




        lpNum++;

    }//end while

    LLast=LCurr;
    y0Last=y0Curr;
    z0Last=z0Curr;
    y1Last=y1Curr;


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
/// @param y0Init
/// @param z0Init
/// @param y1Init
/// @param U_ptr
/// @param L_ptr
/// @param y0_ptr
/// @param z0_ptr
/// @param y1_ptr
void version1DistRelative::executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const double &LInit,
                                              const double &y0Init, const double &z0Init, const double &y1Init,double * U_ptr,double *L_ptr,double *y0_ptr, double *z0_ptr, double *y1_ptr) {


    if (dataNumTotal <= dataNumInEq) {
        return;
    }

    unsigned long long remainingDataNum = this->dataNumTotal - this->dataNumInEq;
    unsigned long long remainingLoopNum = remainingDataNum * lag;
    std::cout << "remainingDataNum=" << remainingDataNum << std::endl;

    std::cout << "remainingLoopNum=" << remainingLoopNum << std::endl;
    unsigned long long lastLoopNum = remainingLoopNum % writeInterval;

    std::cout<<"lastLoopNum="<<lastLoopNum<<std::endl;




    double LCurr=LInit;
    double y0Curr=y0Init;
    double z0Curr=z0Init;
    double y1Curr=y1Init;



    double UCurr = (*potFuncPtr)(LCurr,y0Curr,z0Curr,y1Curr);
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

    std::string out_y0_AllBinSubDir=outDir+"y0_AllBin/";
    std::string out_z0_AllBinSubDir=outDir+"z0_AllBin/";

    std::string out_y1_AllBinSubDir=outDir+"y1_AllBin/";



    const auto tMCStart{std::chrono::steady_clock::now()};
//remainingLoopNum-lastLoopNum
    unsigned long long lpNum = 0;
    while (lpNum < remainingLoopNum - lastLoopNum) {


        //propose a move
        double LNext;
        double y0Next;
        double z0Next;
        double y1Next;

        proposal(LCurr,y0Curr,z0Curr,y1Curr,LNext,y0Next,z0Next,y1Next);
        double r = acceptanceRatio(LCurr,y0Curr,z0Curr,y1Curr,LNext,y0Next,z0Next,y1Next);
        double u = distUnif01(e2);
        if (u <= r) {
            LCurr=LNext;
            y0Curr=y0Next;
            z0Curr=z0Next;
            y1Curr=y1Next;

            UCurr = (*potFuncPtr)(LCurr,y0Curr,z0Curr,y1Curr);
        }//end of accept-reject

        unsigned long long indOfArrayU = lpNum % writeInterval;


        U_ptr[indOfArrayU] = UCurr;
        L_ptr[indOfArrayU]=LCurr;
        y0_ptr[indOfArrayU]=y0Curr;
        z0_ptr[indOfArrayU]=z0Curr;
        y1_ptr[indOfArrayU]=y1Curr;





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

            std::string out_y0_BinFileName=out_y0_AllBinSubDir+filenameMiddle+".y0_All.bin";
            save_to_bin_file(y0_ptr,sizeOfArrayU,out_y0_BinFileName);

            std::string out_z0_BinFileName=out_z0_AllBinSubDir+filenameMiddle+".z0_All.bin";
            save_to_bin_file(z0_ptr,sizeOfArrayU,out_z0_BinFileName);

            std::string out_y1_BinFileName=out_y1_AllBinSubDir+filenameMiddle+".y1_All.bin";
            save_to_bin_file(y1_ptr,sizeOfArrayU,out_y1_BinFileName);



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
        double y0Next;
        double z0Next;
        double y1Next;

        proposal(LCurr,y0Curr,z0Curr,y1Curr,LNext,y0Next,z0Next,y1Next);
        double r = acceptanceRatio(LCurr,y0Curr,z0Curr,y1Curr,LNext,y0Next,z0Next,y1Next);
        double u = distUnif01(e2);
        if (u <= r) {
            LCurr=LNext;
            y0Curr=y0Next;
            z0Curr=z0Next;
            y1Curr=y1Next;

            UCurr = (*potFuncPtr)(LCurr,y0Curr,z0Curr,y1Curr);
        }//end of accept-reject


        U_ptr[lpFinal]=UCurr;
        L_ptr[lpFinal]=LCurr;
        y0_ptr[lpFinal]=y0Curr;
        z0_ptr[lpFinal]=z0Curr;
        y1_ptr[lpFinal]=y1Curr;





    }//end of final for loops

    std::string filenameMiddle = "loopEnd" + std::to_string(lpNum+lastLoopNum-1);
    std::string outUPicleFileName = outUAllPickleSubDir + filenameMiddle + ".UAll.pkl";
    std::string outUBinFileName = outUAllBinSubDir + filenameMiddle + ".UAll.bin";

    save_array_to_pickle(U_ptr,lastLoopNum,outUPicleFileName);
    save_to_bin_file(U_ptr,lastLoopNum,outUBinFileName);

    std::string out_L_BinFileName=out_L_AllBinSubDir+filenameMiddle+".L.bin";
    save_to_bin_file(L_ptr,lastLoopNum,out_L_BinFileName);

    std::string out_y0_BinFileName=out_y0_AllBinSubDir+filenameMiddle+".y0_All.bin";
    save_to_bin_file(y0_ptr,lastLoopNum,out_y0_BinFileName);

    std::string out_z0_BinFileName=out_z0_AllBinSubDir+filenameMiddle+".z0_All.bin";
    save_to_bin_file(z0_ptr,lastLoopNum,out_z0_BinFileName);

    std::string out_y1_BinFileName=out_y1_AllBinSubDir+filenameMiddle+".y1_All.bin";
    save_to_bin_file(y1_ptr,lastLoopNum,out_y1_BinFileName);



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