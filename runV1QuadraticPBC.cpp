#include "./version1/quadraticPBC/version1QuadraticPBC2Atom.hpp"

//running version 1, Lennard-Jones+ quartic+PBC

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }
    double T = std::stod(argv[1]);
    int rowNum=std::stoi(argv[2]);
    double  a1;
    double a2;
    double c1;
    double c2;


    version1Quadratic::parseCSV(rowNum,a1,a2,c1,c2);
    double r0=1;
//    std::cout<<"r0="<<r0<<std::endl;
//    auto qFunc=quadratic(a1,a2,c1,c2,r0);
    unsigned long long cellNum = 10;

    auto v1Obj=version1Quadratic(rowNum,T,cellNum,std::make_shared<quadratic>(a1,a2,c1,c2,r0));

    unsigned long long lag=0;
    unsigned long long totalLoopEq=0;
    bool eq=false;
    bool same= false;
    arma::dcolvec last_xA;
    arma::dcolvec last_xB;
    double last_L;
    std::shared_ptr<double[]> U_ptr;
    std::shared_ptr<double[]> L_ptr;
    std::shared_ptr<double[]> xA_ptr;
    std::shared_ptr<double[]> xB_ptr;

    try{
        U_ptr=std::shared_ptr<double[]>(new double[version1Quadratic::loopMax], std::default_delete<double[]>());
        L_ptr=std::shared_ptr<double[]>(new double[version1Quadratic::loopMax], std::default_delete<double[]>());;
        xA_ptr=std::shared_ptr<double[]>(new double[version1Quadratic::loopMax*cellNum], std::default_delete<double[]>());
        xB_ptr=std::shared_ptr<double[]>(new double[version1Quadratic::loopMax*cellNum], std::default_delete<double[]>());
    }
    catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    v1Obj.readEqMc(lag,totalLoopEq,eq,same,last_xA,last_xB,last_L,U_ptr.get(),L_ptr.get(),xA_ptr.get(),xB_ptr.get());

    std::cout<<"after reacEqMc: equilibrium="<<eq<<std::endl;

    if(!same and lag>0 and eq){
        v1Obj.executionMCAfterEq(lag,totalLoopEq,last_xA,last_xB,last_L,U_ptr.get(),L_ptr.get(),xA_ptr.get(),xB_ptr.get());
    }

    return 0;

}