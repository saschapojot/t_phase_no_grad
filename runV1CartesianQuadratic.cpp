#include "./version1/cartesian/version1CartesianQuadratic.hpp"

//running version 1, Lennard-Jones+ quartic+PBC

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }
    double T = std::stod(argv[1]);
    int rowNum = std::stoi(argv[2]);
    double a1;
    double a2;
    double c1;
    double c2;
    double mA;
    double mB;


    version1CartesianQuadratic::parseCSV(rowNum, a1, a2, c1, c2, mA, mB);

//    std::cout<<"r0="<]<r0<<std::endl;
//    auto qFunc=quadratic(a1,a2,c1,c2,r0);
    unsigned long long cellNum = 2;

    auto v1Obj = version1CartesianQuadratic(rowNum, T, cellNum, std::make_shared<quadraticCartesian>(a1, a2, c1, c2, mA, mB));

    unsigned long long lag = 0;
    unsigned long long totalLoopEq = 0;
    bool eq = false;
    bool same = false;

    double last_L;
    double last_x0A;
    double last_x0B;
    double last_x1A;
    double last_x1B;
    std::shared_ptr<double[]> U_ptr;
    std::shared_ptr<double[]> L_ptr;
    std::shared_ptr<double[]> x0A_ptr;
    std::shared_ptr<double[]> x0B_ptr;
    std::shared_ptr<double[]> x1A_ptr;
    std::shared_ptr<double[]> x1B_ptr;

    try {
        U_ptr = std::shared_ptr<double[]>(new double[version1CartesianQuadratic::loopMax],
                                          std::default_delete<double[]>());
        L_ptr = std::shared_ptr<double[]>(new double[version1CartesianQuadratic::loopMax],
                                          std::default_delete<double[]>());
        x0A_ptr = std::shared_ptr<double[]>(new double[version1CartesianQuadratic::loopMax],
                                            std::default_delete<double[]>());
        x0B_ptr = std::shared_ptr<double[]>(new double[version1CartesianQuadratic::loopMax],
                                            std::default_delete<double[]>());
        x1A_ptr = std::shared_ptr<double[]>(new double[version1CartesianQuadratic::loopMax],
                                            std::default_delete<double[]>());
        x1B_ptr = std::shared_ptr<double[]>(new double[version1CartesianQuadratic::loopMax],
                                            std::default_delete<double[]>());

    }
    catch (const std::bad_alloc &e) {
        std::cerr << "Memory allocation error: " << e.what() << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    v1Obj.readEqMc(lag, totalLoopEq, eq, same, last_L,last_x0A, last_x0B, last_x1A, last_x1B, U_ptr.get(),
                   L_ptr.get(),x0A_ptr.get(), x0B_ptr.get(), x1A_ptr.get(), x1B_ptr.get());

    std::cout << "after reacEqMc: equilibrium=" << eq << std::endl;

    if (!same and lag > 0 and eq) {
        v1Obj.executionMCAfterEq(lag, totalLoopEq, last_L,last_x0A, last_x0B, last_x1A, last_x1B, U_ptr.get(),
                                 L_ptr.get(),x0A_ptr.get(), x0B_ptr.get(), x1A_ptr.get(), x1B_ptr.get());
    }

    return 0;

}