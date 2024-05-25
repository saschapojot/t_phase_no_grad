#include <iostream>
#include "version1/LJPotPBC/version1LJPotPBC2Atom.hpp"



int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }
    double T = std::stod(argv[1]);
    int rowNum=std::stoi(argv[2]);
    double  alpha1;
    double alpha2;
    double beta1;
    double beta2;

    double p1;
    double q1;

    double p2;
    double  q2;
    double r0=-1;
    version1dLJPot2Atom::parseCSV(rowNum,alpha1,beta1,p1,q1,alpha2,beta2,p2,q2,r0);

//    std::cout<<"r0="<<r0<<std::endl;
    auto LJFunc=LJPotPBC(alpha1,alpha2,beta1,beta2,p1,p2,q1,q2,r0);
    unsigned long long cellNum = 10;
    std::cout.precision(11);
    auto v1Obj=version1dLJPot2Atom(rowNum,T,cellNum,std::make_shared<LJPotPBC>(alpha1,alpha2,beta1,beta2,p1,p2,q1,q2,r0));

}
