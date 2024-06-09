#include "version1/cartesian/parseBinPBC.hpp"



int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }


    int rowNum=std::stoi(argv[1]);
    int whichT=std::stoi(argv[2]);
    const auto tStart{std::chrono::steady_clock::now()};
    unsigned long long cellNum=2;
    auto rd=reader(rowNum,whichT,cellNum);
    rd.searchFiles();
    rd.sortFiles();
    rd.parseSummary();
    rd.data2json();


    const auto tEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tEnd - tStart};
    std::cout<< "parsing time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;







}