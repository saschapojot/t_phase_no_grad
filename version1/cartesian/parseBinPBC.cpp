//
// Created by polya on 5/26/24.
//

#include "parseBinPBC.hpp"
void reader::searchFiles() {
    this->UPath = this->TDir + "/UAllBin/";

    this->LPath = this->TDir + "/L_AllBin/";
    this->x0APath=this->TDir + "/x0A_AllBin/";
    this->x0BPath = this->TDir + "/x0B_AllBin/";
    this->x1APath = this->TDir + "/x1A_AllBin/";
    this->x1BPath = this->TDir + "/x1B_AllBin/";
    for (const auto &entry: fs::directory_iterator(UPath)) {
        this->UFilesAll.push_back(entry.path().string());
    }

    for(const auto& entry:fs::directory_iterator(LPath)){
        this->LFilesAll.push_back(entry.path().string());
    }
    for(const auto& entry:fs::directory_iterator(x0APath)){
        this->x0AFilesAll.push_back(entry.path().string());
    }


    for (const auto &entry: fs::directory_iterator(x0BPath)) {
        this->x0BFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(x1APath)) {
        this->x1AFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(x1BPath)) {
        this->x1BFilesAll.push_back(entry.path().string());
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
    this->sorted_LFilesAll.push_back(this->LPath+"/loopStart0ReachEq.L.bin");
    std::vector<std::string> LEndSorted=this->sortOneDir(this->LFilesAll);
    this->sorted_LFilesAll.insert(this->sorted_LFilesAll.end(),LEndSorted.begin(),LEndSorted.end());
//    printVec(sorted_xAFilesAll);

    this->sorted_x0AFilesAll.push_back(this->x0APath+"/loopStart0ReachEq.x0A_All.bin");
    std::vector<std::string> x0AEndSorted=this->sortOneDir(this->x0AFilesAll);
    this->sorted_x0AFilesAll.insert(this->sorted_x0AFilesAll.end(),x0AEndSorted.begin(),x0AEndSorted.end());

    this->sorted_x0BFilesAll.push_back(this->x0BPath+"/loopStart0ReachEq.x0B_All.bin");
    std::vector<std::string> x0BEndSorted=this->sortOneDir(this->x0BFilesAll);
    this->sorted_x0BFilesAll.insert(this->sorted_x0BFilesAll.end(),x0BEndSorted.begin(),x0BEndSorted.end());
//    printVec(sorted_xBFilesAll);

    this->sorted_x1AFilesAll.push_back(this->x1APath+"/loopStart0ReachEq.x1A_All.bin");
    std::vector<std::string> x1AEndSorted=this->sortOneDir(this->x1AFilesAll);
    this->sorted_x1AFilesAll.insert(this->sorted_x1AFilesAll.end(),x1AEndSorted.begin(),x1AEndSorted.end());

    this->sorted_x1BFilesAll.push_back(this->x1BPath+"/loopStart0ReachEq.x1B_All.bin");
    std::vector<std::string> x1BEndSorted=this->sortOneDir(this->x1BFilesAll);
    this->sorted_x1BFilesAll.insert(this->sorted_x1BFilesAll.end(),x1BEndSorted.begin(),x1BEndSorted.end());
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
            lagEst=static_cast<unsigned long long>(static_cast<double >(lagEst));
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
    unsigned long long UNumMax=version1CartesianQuadratic::loopMax+(UFileNum-1)*version1CartesianQuadratic::loopToWrite;

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


unsigned long long reader::load_L(){

    size_t LFileNum=sorted_LFilesAll.size();

    unsigned long long LNumMax=version1CartesianQuadratic::loopMax+(LFileNum-1)*version1CartesianQuadratic::loopToWrite;

    unsigned long long LSelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(LNumMax))/(static_cast<double>(lagEst))));

    this->L_selected=std::shared_ptr<double[]>(new double[LSelectedNum],
                                              std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & LFileName:this->sorted_LFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(LFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            L_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}


unsigned long long  reader::load_x0A(){
    size_t x0AFileNum=sorted_x0AFilesAll.size();
    unsigned long long x0ANumMax=version1CartesianQuadratic::loopMax+(x0AFileNum-1)*version1CartesianQuadratic::loopToWrite;
    unsigned long long x0ASelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(x0ANumMax))/(static_cast<double>(lagEst))));

    this->x0A_selected=std::shared_ptr<double[]>(new double[x0ASelectedNum],
                                                     std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;

    for(const std::string & x0AFileName:this->sorted_x0AFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(x0AFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            x0A_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }

    return counter;

}

unsigned long long reader::load_x0B(){

    size_t x0BFileNum=sorted_x0BFilesAll.size();
    unsigned long long x0BNumMax=version1CartesianQuadratic::loopMax+(x0BFileNum-1)*version1CartesianQuadratic::loopToWrite;
    unsigned long long x0BSelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(x0BNumMax))/(static_cast<double>(lagEst))));


    this->x0B_selected=std::shared_ptr<double[]>(new double[x0BSelectedNum],
                                               std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & x0BFileName:this->sorted_x0BFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(x0BFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            x0B_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}



unsigned long long reader::load_x1A(){

    size_t x1AFileNum=sorted_x1AFilesAll.size();
    unsigned long long x1ANumMax=version1CartesianQuadratic::loopMax+(x1AFileNum-1)*version1CartesianQuadratic::loopToWrite;

    unsigned long long x1ASelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(x1ANumMax))/(static_cast<double>(lagEst))));

    this->x1A_selected=std::shared_ptr<double[]>(new double[x1ASelectedNum],
                                                     std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & x1AFileName:this->sorted_x1AFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(x1AFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            x1A_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}


unsigned long long reader::load_x1B(){
    size_t x1BFileNum=sorted_x1BFilesAll.size();
    unsigned long long x1BNumMax=version1CartesianQuadratic::loopMax+(x1BFileNum-1)*version1CartesianQuadratic::loopToWrite;
    unsigned long long x1BSelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(x1BNumMax))/(static_cast<double>(lagEst))));

    this->x1B_selected=std::shared_ptr<double[]>(new double[x1BSelectedNum],
                                                     std::default_delete<double[]>());

    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & x1BFileName:this->sorted_x1BFilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(x1BFileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            x1B_selected[counter]=db_inOneFile[i];
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

    //L
    const auto tRead_LStart{std::chrono::steady_clock::now()};
    unsigned long long L_size = this->load_L();
    const auto tRead_LEnd{std::chrono::steady_clock::now()};

    const std::chrono::duration<double> elapsed_LsecondsAll{tRead_LEnd - tRead_LStart};
    std::cout << "L_size=" << L_size << std::endl;
    std::cout << "read L time: " << elapsed_LsecondsAll.count() / 3600.0 << " h" << std::endl;
    //write L
    std::string LJsonPath = jsonPath + "/jsonL/";
    if (!fs::is_directory(LJsonPath) || !fs::exists(LJsonPath)) {
        fs::create_directories(LJsonPath);
    }
    std::string LJsonFile = LJsonPath + "/LData.json";
    boost::json::object obj_L;
    boost::json::array arr_L;
    for (unsigned long long i = 0; i < L_size; i++) {
        arr_L.push_back(L_selected[i]);

    }
    obj_L["L"] = arr_L;
    std::ofstream ofsL(LJsonFile);
    std::string LStr = boost::json::serialize(obj_L);
    ofsL << LStr << std::endl;
    ofsL.close();

    //x

    //x0A
    const auto tRead_x0AStart{std::chrono::steady_clock::now()};
    unsigned long long x0A_size = this->load_x0A();
    const auto tRead_x0AEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_x0AsecondsAll{tRead_x0AEnd - tRead_x0AStart};
    std::cout << "x0A_size=" << x0A_size << std::endl;
    std::cout << "read x0A time: " << elapsed_x0AsecondsAll.count() / 3600.0 << " h" << std::endl;

    //x0B
    const auto tRead_x0BStart{std::chrono::steady_clock::now()};
    unsigned long long x0B_size = this->load_x0B();
    const auto tRead_x0BEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_x0BsecondsAll{tRead_x0BEnd - tRead_x0BStart};
    std::cout << "x0B_size=" << x0B_size << std::endl;
    std::cout << "read x0B time: " << elapsed_x0BsecondsAll.count() / 3600.0 << " h" << std::endl;


    //x1A
    const auto tRead_x1AStart{std::chrono::steady_clock::now()};
    unsigned long long x1A_size = this->load_x1A();
    const auto tRead_x1AEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_x1AsecondsAll{tRead_x1AEnd - tRead_x1AStart};
    std::cout << "x1A_size=" << x1A_size << std::endl;
    std::cout << "read x1A time: " << elapsed_x1AsecondsAll.count() / 3600.0 << " h" << std::endl;

    //x1B

    const auto tRead_x1BStart{std::chrono::steady_clock::now()};
    unsigned long long x1B_size = this->load_x1B();
    const auto tRead_x1BEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_x1BsecondsAll{tRead_x1BEnd - tRead_x1BStart};
    std::cout << "x1B_size=" << x1B_size << std::endl;
    std::cout << "read x1B time: " << elapsed_x1BsecondsAll.count() / 3600.0 << " h" << std::endl;

    //write x
    std::string xJsonPath = jsonPath + "/jsonx/";
    if (!fs::is_directory(xJsonPath) || !fs::exists(xJsonPath)) {
        fs::create_directories(xJsonPath);
    }
    std::string xJsonFile = xJsonPath + "/xData.json";

    boost::json::object obj_x;

    boost::json::array arr_x0A;
    boost::json::array arr_x0B;
    boost::json::array arr_x1A;
    boost::json::array arr_x1B;

    for (unsigned long long i = 0; i < x0A_size; i++) {
        arr_x0A.push_back(x0A_selected[i]);

    }

    for (unsigned long long i = 0; i < x0B_size; i++) {
        arr_x0B.push_back(x0B_selected[i]);

    }

    for (unsigned long long i = 0; i < x1A_size; i++) {
        arr_x1A.push_back(x1A_selected[i]);

    }

    for (unsigned long long i = 0; i < x1B_size; i++) {
        arr_x1B.push_back(x1B_selected[i]);

    }
    obj_x["x0A"] = arr_x0A;
    obj_x["x0B"] = arr_x0B;
    obj_x["x1A"] = arr_x1A;
    obj_x["x1B"] = arr_x1B;
    std::ofstream ofsx(xJsonFile);
    std::string xStr = boost::json::serialize(obj_x);
    ofsx << xStr << std::endl;
    ofsx.close();





}



