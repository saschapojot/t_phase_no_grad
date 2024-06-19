#include "./version1/distRelative/version1DistRelative.hpp"

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


    version1DistRelative::parseCSV(rowNum, a1, a2, c1, c2, mA, mB);

//    std::cout<<"r0="<]<r0<<std::endl;
//    auto qFunc=quadratic(a1,a2,c1,c2,r0);
    unsigned long long cellNum = 2;

    auto v1Obj = version1DistRelative(rowNum, T, cellNum, std::make_shared<quadraticDistRelative>(a1, a2, c1, c2, mA, mB));

    unsigned long long lag = 0;
    unsigned long long totalLoopEq = 0;
    bool eq = false;
    bool same = false;

    double last_L;
    double last_y0;
    double last_z0;
    double last_y1;

    std::shared_ptr<double[]> U_ptr;
    std::shared_ptr<double[]> L_ptr;
    std::shared_ptr<double[]> y0_ptr;
    std::shared_ptr<double[]> z0_ptr;
    std::shared_ptr<double[]> y1_ptr;


    try {
        U_ptr = std::shared_ptr<double[]>(new double[version1DistRelative::loopMax],
                                          std::default_delete<double[]>());
        L_ptr = std::shared_ptr<double[]>(new double[version1DistRelative::loopMax],
                                          std::default_delete<double[]>());
        y0_ptr = std::shared_ptr<double[]>(new double[version1DistRelative::loopMax],
                                            std::default_delete<double[]>());
        z0_ptr = std::shared_ptr<double[]>(new double[version1DistRelative::loopMax],
                                            std::default_delete<double[]>());
        y1_ptr = std::shared_ptr<double[]>(new double[version1DistRelative::loopMax],
                                            std::default_delete<double[]>());


    }
    catch (const std::bad_alloc &e) {
        std::cerr << "Memory allocation error: " << e.what() << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    v1Obj.readEqMc(lag, totalLoopEq, eq, same, last_L,last_y0, last_z0, last_y1, U_ptr.get(),
                   L_ptr.get(),y0_ptr.get(), z0_ptr.get(), y1_ptr.get());

    std::cout << "after reacEqMc: equilibrium=" << eq << std::endl;

    if (!same and lag > 0 and eq) {
        v1Obj.executionMCAfterEq(lag, totalLoopEq, last_L,last_y0, last_z0, last_y1, U_ptr.get(),
                                 L_ptr.get(),y0_ptr.get(), z0_ptr.get(), y1_ptr.get());
    }

    return 0;

}