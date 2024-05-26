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
    this->sorted_UFilesAll.push_back(this->UPath+"/loopStart0ReachEqUAll.bin");
    std::vector<std::string> UEndSorted=this->sortOneDir(this->UFilesAll);
    this->sorted_UFilesAll.insert(this->sorted_UFilesAll.end(),UEndSorted.begin(),UEndSorted.end());

    this->sorted_xAFilesAll.push_back(this->xAPath+"/loopStart0ReachEq.xA_All.bin");
    std::vector<std::string> xAEndSorted=this->sortOneDir(this->xAFilesAll);
    this->sorted_xAFilesAll.insert(this->sorted_xAFilesAll.end(),xAEndSorted.begin(),xAEndSorted.end());

    this->sorted_xBFilesAll.push_back(this->xBPath+"/loopStart0ReachEq.xB_All.bin");
    std::vector<std::string> xBEndSorted=this->sortOneDir(this->xBFilesAll);
    this->sorted_xBFilesAll.insert(this->sorted_xBFilesAll.end(),xBEndSorted.begin(),xBEndSorted.end());



}

void reader::parseSummary(){





}