//
// Created by polya on 5/26/24.
//

#include "parseBinPBC.hpp"
///UAll, xA_All, xB_All folder's files
void reader::searchFiles() {
    this->UPath = this->TDir + "/UAllBin/";

    this->xAPath = this->TDir + "/xA_AllBin/";
    this->xBPath = this->TDir + "/xB_AllBin/";
    for (const auto &entry: fs::directory_iterator(UPath)) {
        this->UFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(xAPath)) {
        this->xAFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(xBPath)) {
        this->xBFilesAll.push_back(entry.path().string());
    }


}

///
/// @param path the path containing xml files
/// @return sorted bin files by end loop
std::vector<std::string> reader::sortOneDir(const std::vector<std::string> &allFiles) {
    std::vector<unsigned long long> endLoopsAll;
    std::vector<std::string> filesNamedEnd;
    for (const std::string &name: allFiles) {
        std::regex endPattern("loopEnd(\\d+)");
        std::smatch matchPattern;
        if (std::regex_search(name, matchPattern, endPattern)) {
            filesNamedEnd.push_back(name);
            endLoopsAll.push_back(std::stoull(matchPattern.str(1)));
        }
    }

    std::vector<size_t> inds=this->argsort<unsigned long long>(endLoopsAll);

    std::vector<std::string> sortedFiles;
    for(const auto&i:inds){
        sortedFiles.push_back(filesNamedEnd[i]);
    }
    return sortedFiles;

}

///sort files by starting loop
void reader::sortFiles(){
    this->sorted_UFilesAll.push_back(this->UPath+"/loopStart0ReachEq.UAll.bin");
    std::vector<std::string> UEndSorted=this->sortOneDir(this->UFilesAll);
    this->sorted_UFilesAll.insert(this->sorted_UFilesAll.end(),UEndSorted.begin(),UEndSorted.end());
//    printVec(sorted_UFilesAll);
    this->sorted_xAFilesAll.push_back(this->xAPath+"/loopStart0ReachEq.xA_All.bin");
    std::vector<std::string> xAEndSorted=this->sortOneDir(this->xAFilesAll);
    this->sorted_xAFilesAll.insert(this->sorted_xAFilesAll.end(),xAEndSorted.begin(),xAEndSorted.end());
//    printVec(sorted_xAFilesAll);
    this->sorted_xBFilesAll.push_back(this->xBPath+"/loopStart0ReachEq.xB_All.bin");
    std::vector<std::string> xBEndSorted=this->sortOneDir(this->xBFilesAll);
    this->sorted_xBFilesAll.insert(this->sorted_xBFilesAll.end(),xBEndSorted.begin(),xBEndSorted.end());
//    printVec(sorted_xBFilesAll);


}

void reader::parseSummary(){


    std::string smrPath = TDir + "/summary.txt";
    std::regex lagPattern("lag=([+-]?\\d+)");
    std::regex ctStartPattern("nEqCounterStart=([+-]?\\d+)");

    std::smatch matchLag;
    std::smatch matchCtStart;

    std::ifstream smrIn(smrPath);
    for (std::string line; std::getline(smrIn, line);) {

        //extract lag value
        if (std::regex_search(line, matchLag, lagPattern)) {
            this->lagEst = std::stoull(matchLag.str(1));
            lagEst=static_cast<unsigned long long>(static_cast<double >(lagEst)*1.5);
            std::cout << "lagEst=" << lagEst << std::endl;
        }

        if(std::regex_search(line,matchCtStart,ctStartPattern)){
            this->nEqCounterStart=std::stoull(matchCtStart.str(1));
            std::cout<<"nEqCounterStart="<<nEqCounterStart<<std::endl;
        }

    }//end readline for


}


///
/// @param filename file name
/// @param values values in file
/// @param number_of_values number of values
bool reader::loadMsgFile(const std::string& filename, std::shared_ptr<double[]>& values, size_t& number_of_values){
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    // Read the file content into a buffer
    std::vector<char> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    // Unpack the data to a MessagePack object
    msgpack::object_handle oh = msgpack::unpack(buffer.data(), buffer.size());
    msgpack::object obj = oh.get();

    // Check if the data is an array
    if (obj.type != msgpack::type::ARRAY) {
        std::cerr << "Data is not an array" << std::endl;
        return false;
    }

    // Get the number of values
    number_of_values = obj.via.array.size;

    // Allocate memory for the values
//    values = std::shared_ptr<double[]>(new double[number_of_values]);

    // Load the values into the preallocated shared_ptr<double[]>
    for (size_t i = 0; i < number_of_values; ++i) {
//        if (obj.via.array.ptr[i].type != msgpack::type::FLOAT64) {
//            std::cerr << "Element " << i << " type: " << obj.via.array.ptr[i].type << " is not a double value" << std::endl;
//            std::cerr<<"Element "<< i<<" = "<<obj.via.array.ptr[i]<<std::endl;
//            return false;
//        }
        values[i] = obj.via.array.ptr[i].via.f64;
    }
    return true;

}


unsigned long long reader::loadU(){
    size_t UFileNum=sorted_UFilesAll.size();
    unsigned long long UNumMax=version1Quadratic::loopMax+(UFileNum-1)*version1Quadratic::loopToWrite;

    unsigned long long USelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(UNumMax))/(static_cast<double>(lagEst))));

    this->USelected=std::shared_ptr<double[]>(new double[USelectedNum],
                                              std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & UFileName:this->sorted_UFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(UFileName,UInOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            USelected[counter]=UInOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;



}
unsigned long long reader::load_xA(){
    size_t xA_fileNum=sorted_xAFilesAll.size();
    unsigned long long xA_chunkMax=version1Quadratic::loopMax+(xA_fileNum-1)*version1Quadratic::loopToWrite;
    unsigned long long xA_selectedChunkNum=static_cast<unsigned long long>(std::ceil(static_cast<double >(xA_chunkMax)/(static_cast<double>(lagEst))));

    this->xA_selected=std::shared_ptr<double[]>(new double[xA_selectedChunkNum*cellNum],std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart*cellNum;
    for(const std::string &xAFileName:this->sorted_xAFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(xAFileName,x_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            for(unsigned j=0;j<cellNum;j++){
                xA_selected[counter]=x_inOneFile[i+j];
                counter++;
            }
            i+=lagEst*cellNum;
        }
        unsigned long long rest=lengthInOneFile-(i-lagEst*cellNum);
        pointerStart=lagEst*cellNum-rest;


    }

    return counter;


}

unsigned long long reader::load_xB(){
    size_t xB_fileNum=sorted_xBFilesAll.size();
    unsigned long long xB_chunkMax=version1Quadratic::loopMax+(xB_fileNum-1)*version1Quadratic::loopToWrite;
    unsigned long long xB_selectedChunkNum=static_cast<unsigned long long>(std::ceil(static_cast<double >(xB_chunkMax)/(static_cast<double>(lagEst))));
    this->xB_selected=std::shared_ptr<double[]>(new double[xB_selectedChunkNum*cellNum],std::default_delete<double[]>());

    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart*cellNum;

    for(const std::string &xBFileName:this->sorted_xBFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(xBFileName,x_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            for(unsigned j=0;j<cellNum;j++){
                xB_selected[counter]=x_inOneFile[i+j];
                counter++;
            }
            i+=lagEst*cellNum;
        }
        unsigned long long rest=lengthInOneFile-(i-lagEst*cellNum);
        pointerStart=lagEst*cellNum-rest;
    }
    return counter;

}


///data to json, json as input to plot
void reader::data2json(){

    const auto tReadUStart{std::chrono::steady_clock::now()};

    unsigned  long long USize=this->loadU();

    const auto tReadUEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tReadUEnd - tReadUStart};
    std::cout<<"USize="<<USize<<std::endl;
    std::cout << "read U time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    std::string jsonPath = this->TDir + "/jsonData/";

    //write U
    std::string UJsonPath = jsonPath + "/jsonU/";
    if (!fs::is_directory(UJsonPath) || !fs::exists(UJsonPath)) {
        fs::create_directories(UJsonPath);
    }
    std::string UJsonFile = UJsonPath + "/UData.json";

    boost::json::object objU;
    boost::json::array arrU;
    size_t countU0=0;
    for(unsigned long long i=0;i<USize;i++) {
        arrU.push_back(USelected[i]);
        if (std::abs(USelected[i]) < 1e-5) {
            countU0++;
        }
    }
    std::cout<<"countU0="<<countU0<<std::endl;

    objU["U"] = arrU;
    std::ofstream ofsU(UJsonFile);
    std::string UStr = boost::json::serialize(objU);
    ofsU << UStr << std::endl;
    ofsU.close();

    const auto tRead_xAStart{std::chrono::steady_clock::now()};
    unsigned long long xA_size=this->load_xA();
    const auto tRead_xAEnd{std::chrono::steady_clock::now()};

    const std::chrono::duration<double> elapsed_xAsecondsAll{tRead_xAEnd - tRead_xAStart};
    std::cout<<"xA_size="<<xA_size<<std::endl;
    std::cout << "read xA time: " << elapsed_xAsecondsAll.count() / 3600.0 << " h" << std::endl;
//    double *rawPtr_A=xA_selected.get();
//    arma::dcolvec vecA(rawPtr_A,xA_size, true, true);
    arma::dcolvec vecA(xA_size);
    size_t countxA0=0;
    for(size_t i=0;i<xA_size;i++){
        vecA(i)=xA_selected[i];
        if (std::abs(xA_selected[i])<1e-5){
            countxA0++;
        }
    }
    std::cout<<"countxA0="<<countxA0<<std::endl;

    this->arma_xA=arma::reshape(vecA,cellNum,USize).t();
    std::cout<<"arma_xA shape=("<<arma_xA.n_rows<<", "<<arma_xA.n_cols<<")"<<std::endl;


    const auto tRead_xBStart{std::chrono::steady_clock::now()};
    unsigned long long xB_size=this->load_xB();
    const auto tRead_xBEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_xBsecondsAll{tRead_xBEnd - tRead_xBStart};
    std::cout<<"xB_size="<<xB_size<<std::endl;
    std::cout << "read xB time: " << elapsed_xBsecondsAll.count() / 3600.0 << " h" << std::endl;

//    double *rawPtr_B=xB_selected.get();
//    arma::dcolvec vecB(rawPtr_B,xB_size, true, true);
    size_t countxB0=0;

    arma::dcolvec vecB(xB_size);
    for(size_t i=0;i<xB_size;i++){
        vecB(i)=xB_selected[i];
        if (std::abs(xB_selected[i])<1e-5){
            countxB0++;
        }
    }
    std::cout<<"countxB0="<<countxB0<<std::endl;

    this->arma_xB=arma::reshape(vecB,cellNum,USize).t();
    std::cout<<"arma_xB shape=("<<arma_xB.n_rows<<", "<<arma_xB.n_cols<<")"<<std::endl;


    //write xA, xB
    for(unsigned long long i=0;i<cellNum;i++){
        std::string cellPathTmp = jsonPath + "jsonUnitCell" + std::to_string(i) + "/";
        if (!fs::is_directory(cellPathTmp) || !fs::exists(cellPathTmp)) {
            fs::create_directories(cellPathTmp);
        }

        boost::json::object obj_xAxB;

        std::string cellJsonFile = cellPathTmp + "xAxBData.json";
        boost::json::array oneCol_xA;
        for(unsigned long long j=i;j<xA_size;j+=cellNum){
            oneCol_xA.push_back(xA_selected[j]);
        }
        boost::json::array oneCol_xB;
        for(unsigned long long j=i;j<xB_size;j+=cellNum){
            oneCol_xB.push_back(xB_selected[j]);
        }

        obj_xAxB["xA"] = oneCol_xA;

        obj_xAxB["xB"] = oneCol_xB;
        std::ofstream ofsxAxB(cellJsonFile);
        std::string xAxBStr = boost::json::serialize(obj_xAxB);
        ofsxAxB << xAxBStr << std::endl;
        ofsxAxB.close();


    }
}

///compute the column means of arma_xA, arma_xB
void reader::colmeans(){

    this->E_xARow=arma::mean(arma_xA,0);
    this->E_xBRow=arma::mean(arma_xB,0);

    this->E_xACol=E_xARow.t();
    this->E_xBCol=E_xBRow.t();

    this->E_xA2=E_xACol*E_xARow;
    this->E_xB2=E_xBCol*E_xBRow;


}

///compute correlation functions GAA
void reader::computeGAA(){
    arma::dmat YA = arma::zeros(cellNum, cellNum);

    int Q = arma_xA.n_rows;

    for (int q = 0; q < Q; q++) {
        arma::drowvec rowTmp = arma_xA.row(q);
        arma::dcolvec colTmp = rowTmp.t();
        YA += colTmp * rowTmp;

    }

    double QDB = static_cast<double>(Q);

    YA /= QDB;

//    YA.print("YA:");

    arma::dmat GAA = YA - E_xA2;

    std::string outGAA=TDir+"/GAA.csv";
    std::ofstream ofs(outGAA);

    printMat(GAA,ofs);
    ofs.close();




}

///compute correlation functions GAB
void reader::computeGAB(){
    arma::dmat YAB = arma::zeros(cellNum, cellNum);
    int Q = arma_xA.n_rows;

    for(int q=0;q<Q;q++){
        arma::drowvec rowATmp = arma_xA.row(q);
        arma::dcolvec colATmp=rowATmp.t();

        arma::drowvec rowBTmp=arma_xB.row(q);

        YAB+=colATmp*rowBTmp;

    }

    double QDB = static_cast<double>(Q);

    YAB/=QDB;

    arma::dmat GAB=YAB-E_xACol*E_xBRow;
    std::string outGAB=TDir+"/GAB.csv";
    std::ofstream ofs(outGAB);

    printMat(GAB,ofs);
    ofs.close();



}

///compute correlation functions GBB
void reader::computeGBB(){

    arma::dmat YB = arma::zeros(cellNum, cellNum);

    int Q = arma_xB.n_rows;

    for(int q=0;q<Q;q++){
        arma::drowvec rowTmp = arma_xB.row(q);
        arma::dcolvec colTmp = rowTmp.t();
        YB += colTmp * rowTmp;


    }
    double QDB = static_cast<double>(Q);
    YB/=QDB;

    arma::dmat GBB=YB-E_xB2;
    std::string outGBB=TDir+"/GBB.csv";
    std::ofstream ofs(outGBB);

    printMat(GBB,ofs);
    ofs.close();
}

///compute variance of position variables
void reader::computeVar(){

    arma::drowvec varxAVec=arma::var(arma_xA,0,0);
    arma::drowvec varxBVec=arma::var(arma_xB,0,0);
    std::string jsonPath = this->TDir + "/jsonData/";
//    varxAVec.print("varxAVec");
//    varxBVec.print("varxBVec");
//    std::cout<<"jsonPath="<<jsonPath<<std::endl;

//    ///xA
    std::string jsonVar_xAPath= jsonPath + "/jsonVar_xA/";
    if (!fs::is_directory(jsonVar_xAPath) || !fs::exists(jsonVar_xAPath)) {
        fs::create_directories(jsonVar_xAPath);
    }

    std::string var_xAJsonFile=jsonVar_xAPath+"/var_xA.json";

    boost::json::object objVar_xA;
    boost::json::array arrVar_xA;
    for(auto i=0;i<varxAVec.n_elem;i++){
        arrVar_xA.push_back(varxAVec(i));
    }

    objVar_xA["var_xA"]=arrVar_xA;

    std::ofstream  ofsVar_xA(var_xAJsonFile);
    std::string var_xAStr=boost::json::serialize(objVar_xA);
    ofsVar_xA<<var_xAStr<<std::endl;
    ofsVar_xA.close();

    ///xB
    std::string jsonVar_xBPath= jsonPath + "/jsonVar_xB/";
    if (!fs::is_directory(jsonVar_xBPath) || !fs::exists(jsonVar_xBPath)) {
        fs::create_directories(jsonVar_xBPath);
    }

    std::string var_xBJsonFile=jsonVar_xBPath+"/var_xB.json";
    boost::json::object objVar_xB;
    boost::json::array arrVar_xB;

    for(auto i=0;i<varxBVec.n_elem;i++){
        arrVar_xB.push_back(varxBVec(i));
    }

    objVar_xB["var_xB"]=arrVar_xB;
    std::ofstream  ofsVar_xB(var_xBJsonFile);
    std::string var_xBStr=boost::json::serialize(objVar_xB);
    ofsVar_xB<<var_xBStr<<std::endl;
    ofsVar_xB.close();

}