#include <iostream>
#include "version1/LJPotPBC/version1LJPotPBC2Atom.hpp"

class base{
public: base(const double &inVal){
    this->x=inVal;
}

public:
    virtual double square()=0;


public:
    double x;

};

class derived: base{
public:derived(const double&inVal): base(inVal){}

public:
    double square()override{
        return x*x;
    }

public:
    double x;

};

int main(int argc, char *argv[]) {
    double x=0.1;
    auto sq=derived(x);
    std::cout<<sq.square()<<std::endl;

}
