//
// Created by polya on 5/26/24.
//

#include "parseBinPBC.hpp"
///UAll, xA_All, xB_All folder's files
void reader::searchFiles() {
    this->UPath = this->TDir + "/UAllBin/";

    this->rPath = this->TDir + "/r_AllBin/";
    this->theta0APath=this->TDir + "/theta0A_AllBin/";
    this->theta0BPath = this->TDir + "/theta0B_AllBin/";
    this->theta1APath = this->TDir + "/theta1A_AllBin/";
    this->theta1BPath = this->TDir + "/theta1B_AllBin/";
    for (const auto &entry: fs::directory_iterator(UPath)) {
        this->UFilesAll.push_back(entry.path().string());
    }

    for(const auto& entry:fs::directory_iterator(rPath)){
        this->rFilesAll.push_back(entry.path().string());
    }
    for(const auto& entry:fs::directory_iterator(theta0APath)){
        this->theta0AFilesAll.push_back(entry.path().string());
    }


    for (const auto &entry: fs::directory_iterator(theta0BPath)) {
        this->theta0BFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(theta1APath)) {
        this->theta1AFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(theta1BPath)) {
        this->theta1BFilesAll.push_back(entry.path().string());
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
    this->sorted_rFilesAll.push_back(this->rPath+"/loopStart0ReachEq.r.bin");
    std::vector<std::string> rEndSorted=this->sortOneDir(this->rFilesAll);
    this->sorted_rFilesAll.insert(this->sorted_rFilesAll.end(),rEndSorted.begin(),rEndSorted.end());
//    printVec(sorted_xAFilesAll);

    this->sorted_theta0AFilesAll.push_back(this->theta0APath+"/loopStart0ReachEq.theta0A_All.bin");
    std::vector<std::string> theta0AEndSorted=this->sortOneDir(this->theta0AFilesAll);
    this->sorted_theta0AFilesAll.insert(this->sorted_theta0AFilesAll.end(),theta0AEndSorted.begin(),theta0AEndSorted.end());

    this->sorted_theta0BFilesAll.push_back(this->theta0BPath+"/loopStart0ReachEq.theta0B_All.bin");
    std::vector<std::string> theta0BEndSorted=this->sortOneDir(this->theta0BFilesAll);
    this->sorted_theta0BFilesAll.insert(this->sorted_theta0BFilesAll.end(),theta0BEndSorted.begin(),theta0BEndSorted.end());
//    printVec(sorted_xBFilesAll);

    this->sorted_theta1AFilesAll.push_back(this->theta1APath+"/loopStart0ReachEq.theta1A_All.bin");
    std::vector<std::string> theta1AEndSorted=this->sortOneDir(this->theta1AFilesAll);
    this->sorted_theta1AFilesAll.insert(this->sorted_theta1AFilesAll.end(),theta1AEndSorted.begin(),theta1AEndSorted.end());

    this->sorted_theta1BFilesAll.push_back(this->theta1BPath+"/loopStart0ReachEq.theta1B_All.bin");
    std::vector<std::string> theta1BEndSorted=this->sortOneDir(this->theta1BFilesAll);
    this->sorted_theta1BFilesAll.insert(this->sorted_theta1BFilesAll.end(),theta1BEndSorted.begin(),theta1BEndSorted.end());
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
    unsigned long long UNumMax=version1InsideLQuadratic::loopMax+(UFileNum-1)*version1InsideLQuadratic::loopToWrite;

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


unsigned long long reader::load_r(){

    size_t rFileNum=sorted_rFilesAll.size();

    unsigned long long rNumMax=version1InsideLQuadratic::loopMax+(rFileNum-1)*version1InsideLQuadratic::loopToWrite;

    unsigned long long rSelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(rNumMax))/(static_cast<double>(lagEst))));

    this->r_selected=std::shared_ptr<double[]>(new double[rSelectedNum],
                                              std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & rFileName:this->sorted_rFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(rFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            r_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}


unsigned long long  reader::load_theta0A(){
    size_t theta0AFileNum=sorted_theta0AFilesAll.size();
    unsigned long long theta0ANumMax=version1InsideLQuadratic::loopMax+(theta0AFileNum-1)*version1InsideLQuadratic::loopToWrite;
    unsigned long long theta0ASelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(theta0ANumMax))/(static_cast<double>(lagEst))));

    this->theta0A_selected=std::shared_ptr<double[]>(new double[theta0ASelectedNum],
                                                     std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;

    for(const std::string & theta0AFileName:this->sorted_theta0AFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(theta0AFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            theta0A_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }

    return counter;

}

unsigned long long reader::load_theta0B(){

    size_t theta0BFileNum=sorted_theta0BFilesAll.size();
    unsigned long long theta0BNumMax=version1InsideLQuadratic::loopMax+(theta0BFileNum-1)*version1InsideLQuadratic::loopToWrite;
    unsigned long long theta0BSelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(theta0BNumMax))/(static_cast<double>(lagEst))));


    this->theta0B_selected=std::shared_ptr<double[]>(new double[theta0BSelectedNum],
                                               std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & theta0BFileName:this->sorted_theta0BFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(theta0BFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            theta0B_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}



unsigned long long reader::load_theta1A(){

    size_t theta1AFileNum=sorted_theta1AFilesAll.size();
    unsigned long long theta1ANumMax=version1InsideLQuadratic::loopMax+(theta1AFileNum-1)*version1InsideLQuadratic::loopToWrite;

    unsigned long long theta1ASelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(theta1ANumMax))/(static_cast<double>(lagEst))));

    this->theta1A_selected=std::shared_ptr<double[]>(new double[theta1ASelectedNum],
                                                     std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & theta1AFileName:this->sorted_theta1AFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(theta1AFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            theta1A_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}


unsigned long long reader::load_theta1B(){
    size_t theta1BFileNum=sorted_theta1BFilesAll.size();
    unsigned long long theta1BNumMax=version1InsideLQuadratic::loopMax+(theta1BFileNum-1)*version1InsideLQuadratic::loopToWrite;
    unsigned long long theta1BSelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(theta1BNumMax))/(static_cast<double>(lagEst))));

    this->theta1B_selected=std::shared_ptr<double[]>(new double[theta1BSelectedNum],
                                                     std::default_delete<double[]>());

    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & theta1BFileName:this->sorted_theta1BFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(theta1BFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            theta1B_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}




///data to json, json as input to plot
void reader::data2json() {
    //U
    const auto tReadUStart{std::chrono::steady_clock::now()};

    unsigned long long USize = this->loadU();

    const auto tReadUEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tReadUEnd - tReadUStart};
    std::cout << "USize=" << USize << std::endl;
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
    size_t countU0 = 0;
    for (unsigned long long i = 0; i < USize; i++) {
        arrU.push_back(USelected[i]);
        if (std::abs(USelected[i]) < 1e-5) {
            countU0++;
        }
    }
    std::cout << "countU0=" << countU0 << std::endl;

    objU["U"] = arrU;
    std::ofstream ofsU(UJsonFile);
    std::string UStr = boost::json::serialize(objU);
    ofsU << UStr << std::endl;
    ofsU.close();

    //r
    const auto tRead_rStart{std::chrono::steady_clock::now()};
    unsigned long long r_size = this->load_r();
    const auto tRead_rEnd{std::chrono::steady_clock::now()};

    const std::chrono::duration<double> elapsed_rsecondsAll{tRead_rEnd - tRead_rStart};
    std::cout << "r_size=" << r_size << std::endl;
    std::cout << "read r time: " << elapsed_rsecondsAll.count() / 3600.0 << " h" << std::endl;
    //write r
    std::string rJsonPath = jsonPath + "/jsonr/";
    if (!fs::is_directory(rJsonPath) || !fs::exists(rJsonPath)) {
        fs::create_directories(rJsonPath);
    }
    std::string rJsonFile = rJsonPath + "/rData.json";
    boost::json::object obj_r;
    boost::json::array arr_r;
    for (unsigned long long i = 0; i < r_size; i++) {
        arr_r.push_back(r_selected[i]);

    }
    obj_r["r"] = arr_r;
    std::ofstream ofsr(rJsonFile);
    std::string rStr = boost::json::serialize(obj_r);
    ofsr << rStr << std::endl;
    ofsr.close();

    //theta

    //theta0A
    const auto tRead_theta0AStart{std::chrono::steady_clock::now()};
    unsigned long long theta0A_size = this->load_theta0A();
    const auto tRead_theta0AEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_theta0AsecondsAll{tRead_theta0AEnd - tRead_theta0AStart};
    std::cout << "theta0A_size=" << theta0A_size << std::endl;
    std::cout << "read theta0A time: " << elapsed_theta0AsecondsAll.count() / 3600.0 << " h" << std::endl;

    //theta0B
    const auto tRead_theta0BStart{std::chrono::steady_clock::now()};
    unsigned long long theta0B_size = this->load_theta0B();
    const auto tRead_theta0BEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_theta0BsecondsAll{tRead_theta0BEnd - tRead_theta0BStart};
    std::cout << "theta0B_size=" << theta0B_size << std::endl;
    std::cout << "read theta0B time: " << elapsed_theta0BsecondsAll.count() / 3600.0 << " h" << std::endl;


    //theta1A
    const auto tRead_theta1AStart{std::chrono::steady_clock::now()};
    unsigned long long theta1A_size = this->load_theta1A();
    const auto tRead_theta1AEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_theta1AsecondsAll{tRead_theta1AEnd - tRead_theta1AStart};
    std::cout << "theta1A_size=" << theta1A_size << std::endl;
    std::cout << "read theta1A time: " << elapsed_theta1AsecondsAll.count() / 3600.0 << " h" << std::endl;

    //theta1B

    const auto tRead_theta1BStart{std::chrono::steady_clock::now()};
    unsigned long long theta1B_size = this->load_theta1B();
    const auto tRead_theta1BEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_theta1BsecondsAll{tRead_theta1BEnd - tRead_theta1BStart};
    std::cout << "theta1B_size=" << theta1B_size << std::endl;
    std::cout << "read theta1B time: " << elapsed_theta1BsecondsAll.count() / 3600.0 << " h" << std::endl;

    //write theta
    std::string thetaJsonPath = jsonPath + "/jsontheta/";
    if (!fs::is_directory(thetaJsonPath) || !fs::exists(thetaJsonPath)) {
        fs::create_directories(thetaJsonPath);
    }
    std::string thetaJsonFile = thetaJsonPath + "/thetaData.json";

    boost::json::object obj_theta;

    boost::json::array arr_theta0A;
    boost::json::array arr_theta0B;
    boost::json::array arr_theta1A;
    boost::json::array arr_theta1B;

    for (unsigned long long i = 0; i < theta0A_size; i++) {
        arr_theta0A.push_back(theta0A_selected[i]);

    }

    for (unsigned long long i = 0; i < theta0B_size; i++) {
        arr_theta0B.push_back(theta0B_selected[i]);

    }

    for (unsigned long long i = 0; i < theta1A_size; i++) {
        arr_theta1A.push_back(theta1A_selected[i]);

    }

    for (unsigned long long i = 0; i < theta1B_size; i++) {
        arr_theta1B.push_back(theta1B_selected[i]);

    }
    obj_theta["theta0A"] = arr_theta0A;
    obj_theta["theta0B"] = arr_theta0B;
    obj_theta["theta1A"] = arr_theta1A;
    obj_theta["theta1B"] = arr_theta1B;
    std::ofstream ofstheta(thetaJsonFile);
    std::string thetaStr = boost::json::serialize(obj_theta);
    ofstheta << thetaStr << std::endl;
    ofstheta.close();





}



