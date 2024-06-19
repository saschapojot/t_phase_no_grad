//
// Created by polya on 5/26/24.
//

#include "parseBinPBC.hpp"
void reader::searchFiles() {
    this->UPath = this->TDir + "/UAllBin/";

    this->LPath = this->TDir + "/L_AllBin/";
    this->y0Path=this->TDir + "/y0_AllBin/";
    this->z0Path = this->TDir + "/z0_AllBin/";
    this->y1Path = this->TDir + "/y1_AllBin/";

    for (const auto &entry: fs::directory_iterator(UPath)) {
        this->UFilesAll.push_back(entry.path().string());
    }

    for(const auto& entry:fs::directory_iterator(LPath)){
        this->LFilesAll.push_back(entry.path().string());
    }
    for(const auto& entry:fs::directory_iterator(y0Path)){
        this->y0FilesAll.push_back(entry.path().string());
    }


    for (const auto &entry: fs::directory_iterator(z0Path)) {
        this->z0FilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(y1Path)) {
        this->y1FilesAll.push_back(entry.path().string());
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

    this->sorted_y0FilesAll.push_back(this->y0Path+"/loopStart0ReachEq.y0_All.bin");
    std::vector<std::string> y0EndSorted=this->sortOneDir(this->y0FilesAll);
    this->sorted_y0FilesAll.insert(this->sorted_y0FilesAll.end(),y0EndSorted.begin(),y0EndSorted.end());

    this->sorted_z0FilesAll.push_back(this->z0Path+"/loopStart0ReachEq.z0_All.bin");
    std::vector<std::string> z0EndSorted=this->sortOneDir(this->z0FilesAll);
    this->sorted_z0FilesAll.insert(this->sorted_z0FilesAll.end(),z0EndSorted.begin(),z0EndSorted.end());
//    printVec(sorted_xBFilesAll);

    this->sorted_y1FilesAll.push_back(this->y1Path+"/loopStart0ReachEq.y1_All.bin");
    std::vector<std::string> y1EndSorted=this->sortOneDir(this->y1FilesAll);
    this->sorted_y1FilesAll.insert(this->sorted_y1FilesAll.end(),y1EndSorted.begin(),y1EndSorted.end());


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
    unsigned long long UNumMax=version1DistRelative::loopMax+(UFileNum-1)*version1DistRelative::loopToWrite;

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

    unsigned long long LNumMax=version1DistRelative::loopMax+(LFileNum-1)*version1DistRelative::loopToWrite;

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


unsigned long long  reader::load_y0(){
    size_t y0FileNum=sorted_y0FilesAll.size();
    unsigned long long y0NumMax=version1DistRelative::loopMax+(y0FileNum-1)*version1DistRelative::loopToWrite;
    unsigned long long y0SelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(y0NumMax))/(static_cast<double>(lagEst))));

    this->y0_selected=std::shared_ptr<double[]>(new double[y0SelectedNum],
                                                     std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;

    for(const std::string & y0FileName:this->sorted_y0FilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(y0FileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            y0_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }

    return counter;

}

unsigned long long reader::load_z0(){

    size_t z0FileNum=sorted_z0FilesAll.size();
    unsigned long long z0NumMax=version1DistRelative::loopMax+(z0FileNum-1)*version1DistRelative::loopToWrite;
    unsigned long long z0SelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(z0NumMax))/(static_cast<double>(lagEst))));


    this->z0_selected=std::shared_ptr<double[]>(new double[z0SelectedNum],
                                               std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & z0FileName:this->sorted_z0FilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(z0FileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            z0_selected[counter]=db_inOneFile[i];
            i+=lagEst;
            counter++;

        }

        unsigned long long rest=lengthInOneFile-(i-lagEst);
        pointerStart=lagEst-rest;

    }
    return counter;
}



unsigned long long reader::load_y1(){

    size_t y1FileNum=sorted_y1FilesAll.size();
    unsigned long long y1NumMax=version1DistRelative::loopMax+(y1FileNum-1)*version1DistRelative::loopToWrite;

    unsigned long long y1SelectedNum=static_cast<unsigned long long>(std::ceil((static_cast<double >(y1NumMax))/(static_cast<double>(lagEst))));

    this->y1_selected=std::shared_ptr<double[]>(new double[y1SelectedNum],
                                                     std::default_delete<double[]>());
    unsigned long long counter=0;
    unsigned long long pointerStart=this->nEqCounterStart;
    for(const std::string & y1FileName:this->sorted_y1FilesAll){
        size_t lengthInOneFile=0;
        this->loadMsgFile(y1FileName,db_inOneFile,lengthInOneFile);
        unsigned long long i=pointerStart;
        while(i<lengthInOneFile){
            y1_selected[counter]=db_inOneFile[i];
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



    //y0
    const auto tRead_y0Start{std::chrono::steady_clock::now()};
    unsigned long long y0_size = this->load_y0();
    const auto tRead_y0End{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_y0secondsAll{tRead_y0End - tRead_y0Start};
    std::cout << "y0_size=" << y0_size << std::endl;
    std::cout << "read y0 time: " << elapsed_y0secondsAll.count() / 3600.0 << " h" << std::endl;

    //z0
    const auto tRead_z0Start{std::chrono::steady_clock::now()};
    unsigned long long z0_size = this->load_z0();
    const auto tRead_z0End{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_z0secondsAll{tRead_z0End - tRead_z0Start};
    std::cout << "z0_size=" << z0_size << std::endl;
    std::cout << "read z0 time: " << elapsed_z0secondsAll.count() / 3600.0 << " h" << std::endl;


    //y1
    const auto tRead_y1Start{std::chrono::steady_clock::now()};
    unsigned long long y1_size = this->load_y1();
    const auto tRead_y1End{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_y1secondsAll{tRead_y1End - tRead_y1Start};
    std::cout << "y1_size=" << y1_size << std::endl;
    std::cout << "read y1 time: " << elapsed_y1secondsAll.count() / 3600.0 << " h" << std::endl;


    //write dist
    std::string distJsonPath = jsonPath + "/jsondist/";
    if (!fs::is_directory(distJsonPath) || !fs::exists(distJsonPath)) {
        fs::create_directories(distJsonPath);
    }
    std::string distJsonFile = distJsonPath + "/distData.json";

    boost::json::object obj_dist;

    boost::json::array arr_y0;
    boost::json::array arr_z0;
    boost::json::array arr_y1;


    for (unsigned long long i = 0; i < y0_size; i++) {
        arr_y0.push_back(y0_selected[i]);

    }

    for (unsigned long long i = 0; i < z0_size; i++) {
        arr_z0.push_back(z0_selected[i]);

    }

    for (unsigned long long i = 0; i < y1_size; i++) {
        arr_y1.push_back(y1_selected[i]);

    }


    obj_dist["y0"] = arr_y0;
    obj_dist["z0"] = arr_z0;
    obj_dist["y1"] = arr_y1;

    std::ofstream ofsdist(distJsonFile);
    std::string distStr = boost::json::serialize(obj_dist);
    ofsdist << distStr << std::endl;
    ofsdist.close();





}



