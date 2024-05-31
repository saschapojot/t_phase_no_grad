//
// Created by polya on 5/23/24.
//
#include "version1QuadraticPBC2Atom.hpp"


void version1Quadratic::parseCSV(const int &rowNum, double &a1, double &a2, double &c1, double &c2) {

    std::string filePath = "./version1Input/1d/quadratic/quadraticParams.csv";
    std::string pyFile = "./version1/quadraticPBC/readCSV.py";
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










//std::cout<<"r0="<<r0<<std::endl;



}

///
/// @param cmd python execution string
/// @return signal from the python
std::string version1Quadratic::execPython(const char *cmd) {
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
double version1Quadratic::f(const arma::dcolvec &xA, const arma::dcolvec &xB, const double & L){
    return this->beta * ((*potFuncPtr)(xA, xB,L));



}


///
/// @param xACurr positions of atom A
/// @param xBCurr positions of atom B
/// @param LCurr total length
/// @param zANext proposed positions of atom A
/// @param zBNext proposed positions of atom B
/// @param LNext proposed value of length
void version1Quadratic::proposal(const arma::dcolvec &xACurr, const arma::dcolvec &xBCurr, const double &LCurr,
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
double version1Quadratic::acceptanceRatio(const arma::dcolvec &xA, const arma::dcolvec &xB,const double& LCurr,
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
void version1Quadratic::initPositionsEquiDistance(arma::dcolvec &xAInit, arma::dcolvec &xBInit, double &LInit){
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
void version1Quadratic::readEqMc(unsigned long long &lag, unsigned long long&loopTotal, bool &equilibrium, bool &same, arma::dcolvec &xALast,
                                   arma::dcolvec &xBLast, double &LLast,
                                   double * U_ptr,double *L_ptr, double *xA_ptr, double *xB_ptr ){

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
    std::string out_xA_AllBinSubDir = outDir + "xA_AllBin/";
    std::string out_xB_AllBinSubDir = outDir + "xB_AllBin/";
    std::string outLAllBinSubDir=outDir+"LAllBin/";
    if (!fs::is_directory(outUAllPickleSubDir) || !fs::exists(outUAllPickleSubDir)) {
        fs::create_directories(outUAllPickleSubDir);
    }
    if (!fs::is_directory(outUAllBinSubDir) || !fs::exists(outUAllBinSubDir)) {
        fs::create_directories(outUAllBinSubDir);
    }
    if (!fs::is_directory(out_xA_AllBinSubDir) || !fs::exists(out_xA_AllBinSubDir)) {
        fs::create_directories(out_xA_AllBinSubDir);
    }
    if (!fs::is_directory(out_xB_AllBinSubDir ) || !fs::exists(out_xB_AllBinSubDir )) {
        fs::create_directories(out_xB_AllBinSubDir );
    }

    if (!fs::is_directory(outLAllBinSubDir ) || !fs::exists(outLAllBinSubDir )) {
        fs::create_directories(outLAllBinSubDir );
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

//    std::unique_ptr<double[]> U_ptr;
//    std::unique_ptr<double[]> L_ptr;
//    std::unique_ptr<double[]> xA_ptr;
//    std::unique_ptr<double[]> xB_ptr;
//
//    try{
//        U_ptr=std::make_unique<double[]>(loopMax);
//         L_ptr=std::make_unique<double[]>(loopMax);
//        xA_ptr=std::make_unique<double[]>(loopMax*N);
//         xB_ptr=std::make_unique<double[]>(loopMax*N);
//    }
//    catch (const std::bad_alloc& e) {
//        std::cerr << "Memory allocation error: " << e.what() << std::endl;
//    } catch (const std::exception& e) {
//        std::cerr << "Exception: " << e.what() << std::endl;
//    }


    arma::dcolvec xACurr(N);
    arma::dcolvec xBCurr(N);
    double LCurr=0;

    this->initPositionsEquiDistance(xACurr, xBCurr,LCurr);
    double UCurr = (*potFuncPtr)(xACurr, xBCurr,LCurr);
    int lpNum=0;
    bool active = true;
    const auto tMCStart{std::chrono::steady_clock::now()};



    while(lpNum<this->loopMax and active==true){
        //propose a move
        arma::dcolvec xANext = arma::dcolvec(N);
        arma::dcolvec xBNext = arma::dcolvec(N);
        double LNext;
        proposal(xACurr, xBCurr,LCurr, xANext, xBNext,LNext);
        double r = acceptanceRatio(xACurr, xBCurr,LCurr, xANext, xBNext,LNext);
        double u = distUnif01(e2);

        if (u <= r) {
//                std::cout<<"UCurr="<<UCurr<<std::endl;
            xACurr = xANext;
            xBCurr = xBNext;
            LCurr=LNext;
            UCurr = (*potFuncPtr)(xACurr, xBCurr,LCurr);
//                std::cout<<"UNext="<<UCurr<<std::endl;
        }//end of accept-reject

        U_ptr[lpNum]=UCurr;
        L_ptr[lpNum]=LCurr;
        for(unsigned long  long i=0;i<N;i++){
            xA_ptr[i+lpNum*N]=xACurr(i);
            xB_ptr[i+lpNum*N]=xBCurr(i);

        }//end writing to array

        //write to file every loopToWrite loops, and inquire equilibrium
        if ((lpNum+1)%loopToWrite==0 and lpNum>1){
            unsigned long long sizeOfArray=lpNum+1;
            unsigned long long sizeOfCoords=sizeOfArray*N;
//            unsigned long long lpStart=0;

            std::string filenameMiddle="loopStart0ReachEq";
            std::string outUPickleFileName=outUAllPickleSubDir+filenameMiddle+ ".UAll.pkl";
            std::string outUBinFileName=outUAllBinSubDir+filenameMiddle+".UAll.bin";
            save_array_to_pickle(U_ptr,sizeOfArray,outUPickleFileName);
            save_to_bin_file(U_ptr,sizeOfArray,outUBinFileName);


            std::string out_xA_BinFileName=out_xA_AllBinSubDir+filenameMiddle+".xA_All.bin";
            save_to_bin_file(xA_ptr,sizeOfCoords,out_xA_BinFileName);

            std::string out_xB_BinFileName=out_xB_AllBinSubDir+filenameMiddle+".xB_All.bin";
            save_to_bin_file(xB_ptr,sizeOfCoords,out_xB_BinFileName);

            std::string outLBinFileName=outLAllBinSubDir+filenameMiddle+".LAll.bin";
            save_to_bin_file(L_ptr,sizeOfArray,outLBinFileName);


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

    xALast=xACurr;
    xBLast=xBCurr;
    LLast=LCurr;
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
/// @param lag decorrelation length
/// @param loopEq total loop numbers in reaching equilibrium
/// @param xA_init xA from readEqMc
/// @param xB_init xB from readEqMc
/// @param LInit L from readEqMc
void version1Quadratic::executionMCAfterEq(const unsigned long long &lag, const unsigned long long &loopEq, const arma::dcolvec &xA_init,
                        const arma::dcolvec &xB_init, const double &LInit,
                                             double * U_ptr,double *L_ptr, double *xA_ptr, double *xB_ptr ) {


    if (dataNumTotal <= dataNumInEq) {
        return;
    }

    unsigned long long remainingDataNum = this->dataNumTotal - this->dataNumInEq;
    unsigned long long remainingLoopNum = remainingDataNum * lag;
    std::cout << "remainingDataNum=" << remainingDataNum << std::endl;

    std::cout << "remainingLoopNum=" << remainingLoopNum << std::endl;
    unsigned long long lastLoopNum = remainingLoopNum % writeInterval;

    std::cout<<"lastLoopNum="<<lastLoopNum<<std::endl;


//    std::unique_ptr<double[]> U_ptr;
//    std::unique_ptr<double[]> L_ptr;
//    std::unique_ptr<double[]> xA_ptr;
//    std::unique_ptr<double[]> xB_ptr;
//
//    try {
//        U_ptr = std::make_unique<double[]>(writeInterval);
//        L_ptr = std::make_unique<double[]>(writeInterval);
//        xA_ptr = std::make_unique<double[]>(writeInterval * N);
//        xB_ptr = std::make_unique<double[]>(writeInterval * N);
//    }
//    catch (const std::bad_alloc &e) {
//        std::cerr << "Memory allocation error: " << e.what() << std::endl;
//    } catch (const std::exception &e) {
//        std::cerr << "Exception: " << e.what() << std::endl;
//    }

    arma::dcolvec xACurr(xA_init);
    arma::dcolvec xBCurr(xB_init);

    double LCurr = LInit;
    double UCurr = (*potFuncPtr)(xACurr, xBCurr, LCurr);
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
    std::string out_xA_AllBinSubDir = outDir + "xA_AllBin/";
    std::string out_xB_AllBinSubDir = outDir + "xB_AllBin/";
    std::string outLAllBinSubDir = outDir + "LAllBin/";
    const auto tMCStart{std::chrono::steady_clock::now()};
//remainingLoopNum-lastLoopNum
    unsigned long long lpNum = 0;
    while (lpNum < remainingLoopNum - lastLoopNum) {


        //propose a move
        arma::dcolvec xANext = arma::dcolvec(N);
        arma::dcolvec xBNext = arma::dcolvec(N);
        double LNext = 0;
        proposal(xACurr, xBCurr, LCurr, xANext, xBNext, LNext);
        double r = acceptanceRatio(xACurr, xBCurr, LCurr, xANext, xBNext, LNext);
        double u = distUnif01(e2);
        if (u <= r) {
            xACurr = xANext;
            xBCurr = xBNext;
            LCurr = LNext;
            UCurr = (*potFuncPtr)(xACurr, xBCurr, LCurr);
        }//end of accept-reject

        unsigned long long indOfArrayU = lpNum % writeInterval;


        U_ptr[indOfArrayU] = UCurr;
        L_ptr[indOfArrayU] = LCurr;

        for (unsigned long long i = 0; i < N; i++) {
            xA_ptr[i + indOfArrayU * N] = xACurr(i);
            xB_ptr[i + indOfArrayU * N] = xBCurr(i);

        }//end writing to array


        //write to file every writeInterval loops

        if ((lpNum + 1) % loopToWrite == 0 and lpNum > 1) {
            unsigned long long sizeOfArrayU = indOfArrayU + 1;
            unsigned long long sizeOfCoords = N * sizeOfArrayU;
            std::string filenameMiddle = "loopEnd" + std::to_string(lpNum);
            std::string outUPicleFileName = outUAllPickleSubDir + filenameMiddle + ".UAll.pkl";
            std::string outUBinFileName = outUAllBinSubDir + filenameMiddle + "UAll.bin";
            save_array_to_pickle(U_ptr, sizeOfArrayU, outUPicleFileName);
            save_to_bin_file(U_ptr, sizeOfArrayU, outUBinFileName);

            std::string out_xA_BinFileName = out_xA_AllBinSubDir + filenameMiddle + ".xA_All.bin";
            save_to_bin_file(xA_ptr, sizeOfCoords, out_xA_BinFileName);

            std::string out_xB_BinFileName = out_xB_AllBinSubDir + filenameMiddle + ".xB_All.bin";
            save_to_bin_file(xB_ptr, sizeOfCoords, out_xB_BinFileName);

            std::string outLBinFileName = outLAllBinSubDir + filenameMiddle + ".LAll.bin";
            save_to_bin_file(L_ptr, sizeOfArrayU, outLBinFileName);

            const auto tWriteEnd{std::chrono::steady_clock::now()};

            const std::chrono::duration<double> elapsed_seconds{tWriteEnd - tMCStart};
            std::cout << "loop " << lpNum << std::endl;
            std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;


        }//end write (if)

        lpNum++;
    }//end while

    for(unsigned long long lpFinal=0;lpFinal<lastLoopNum;lpFinal++){

        //propose a move
        arma::dcolvec xANext = arma::dcolvec(N);
        arma::dcolvec xBNext = arma::dcolvec(N);
        double LNext = 0;
        proposal(xACurr, xBCurr, LCurr, xANext, xBNext, LNext);
        double r = acceptanceRatio(xACurr, xBCurr, LCurr, xANext, xBNext, LNext);
        double u = distUnif01(e2);
        if (u <= r) {
            xACurr = xANext;
            xBCurr = xBNext;
            LCurr = LNext;
            UCurr = (*potFuncPtr)(xACurr, xBCurr, LCurr);
        }//end of accept-reject


        U_ptr[lpFinal]=UCurr;
        L_ptr[lpFinal]=LCurr;

        for(unsigned long long i=0;i<N;i++){
            xA_ptr[i+lpFinal*N]=xACurr(i);
            xB_ptr[i+lpFinal*N]=xBCurr(i);

        }//end writing to array


    }//end of final for loops

    std::string filenameMiddle = "loopEnd" + std::to_string(lpNum+lastLoopNum-1);
    std::string outUPicleFileName = outUAllPickleSubDir + filenameMiddle + ".UAll.pkl";
    std::string outUBinFileName = outUAllBinSubDir + filenameMiddle + "UAll.bin";

    save_array_to_pickle(U_ptr,lastLoopNum,outUPicleFileName);
    save_to_bin_file(U_ptr,lastLoopNum,outUBinFileName);

    std::string out_xA_BinFileName = out_xA_AllBinSubDir + filenameMiddle + ".xA_All.bin";
    save_to_bin_file(xA_ptr,lastLoopNum*N,out_xA_BinFileName);
    std::string out_xB_BinFileName = out_xB_AllBinSubDir + filenameMiddle + ".xB_All.bin";

    save_to_bin_file(xB_ptr,lastLoopNum*N,out_xB_BinFileName);

    std::string outLBinFileName = outLAllBinSubDir + filenameMiddle + ".LAll.bin";

    save_to_bin_file(L_ptr,lastLoopNum,outLBinFileName);

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