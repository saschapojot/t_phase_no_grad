#include "version1/quadraticPBC/parseBinPBC.hpp"



int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }


    int rowNum=std::stoi(argv[1]);
    int whichT=std::stoi(argv[2]);
    const auto tStart{std::chrono::steady_clock::now()};
    unsigned long long cellNum=10;
    auto rd=reader(rowNum,whichT,cellNum);
    rd.searchFiles();
    rd.sortFiles();
    rd.parseSummary();
    rd.data2json();
    rd.colmeans();
    rd.computeGAA();
    rd.computeGAB();
    rd.computeGBB();







}