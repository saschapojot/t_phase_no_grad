#include "./version1/insideLQuadratic/version1InsideLQuadratic.hpp"

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


    version1InsideLQuadratic::parseCSV(rowNum, a1, a2, c1, c2, mA, mB);

//    std::cout<<"r0="<]<r0<<std::endl;
//    auto qFunc=quadratic(a1,a2,c1,c2,r0);
    unsigned long long cellNum = 2;

    auto v1Obj = version1InsideLQuadratic(rowNum, T, cellNum, std::make_shared<quadraticInsideL>(a1, a2, c1, c2, mA, mB));


    unsigned long long lag = 0;
    unsigned long long totalLoopEq = 0;
    bool eq = false;
    bool same = false;

    double last_r;
    double last_theta0A;
    double last_theta0B;
    double last_theta1A;
    double last_theta1B;
    std::shared_ptr<double[]> U_ptr;
    std::shared_ptr<double[]> r_ptr;
    std::shared_ptr<double[]> theta0A_ptr;
    std::shared_ptr<double[]> theta0B_ptr;
    std::shared_ptr<double[]> theta1A_ptr;
    std::shared_ptr<double[]> theta1B_ptr;

    try {
        U_ptr = std::shared_ptr<double[]>(new double[version1InsideLQuadratic::loopMax], std::default_delete<double[]>());
        r_ptr = std::shared_ptr<double[]>(new double[version1InsideLQuadratic::loopMax], std::default_delete<double[]>());
        theta0A_ptr = std::shared_ptr<double[]>(new double[version1InsideLQuadratic::loopMax],
                                                std::default_delete<double[]>());
        theta0B_ptr = std::shared_ptr<double[]>(new double[version1InsideLQuadratic::loopMax],
                                                std::default_delete<double[]>());
        theta1A_ptr = std::shared_ptr<double[]>(new double[version1InsideLQuadratic::loopMax],
                                                std::default_delete<double[]>());
        theta1B_ptr = std::shared_ptr<double[]>(new double[version1InsideLQuadratic::loopMax],
                                                std::default_delete<double[]>());

    }
    catch (const std::bad_alloc &e) {
        std::cerr << "Memory allocation error: " << e.what() << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    v1Obj.readEqMc(lag, totalLoopEq, eq, same, last_r,last_theta0A, last_theta0B, last_theta1A, last_theta1B, U_ptr.get(),
                   r_ptr.get(),theta0A_ptr.get(), theta0B_ptr.get(), theta1A_ptr.get(), theta1B_ptr.get());

    std::cout << "after reacEqMc: equilibrium=" << eq << std::endl;

    if (!same and lag > 0 and eq) {
        v1Obj.executionMCAfterEq(lag, totalLoopEq, last_r,last_theta0A, last_theta0B, last_theta1A, last_theta1B, U_ptr.get(),
                                 r_ptr.get(),theta0A_ptr.get(), theta0B_ptr.get(), theta1A_ptr.get(), theta1B_ptr.get());
    }

    return 0;

}