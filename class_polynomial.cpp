#ifndef _POLY_
#define _POLY_
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include "fastSort_Vector.cpp"
#endif

using namespace std;
const int init = 0.0;

class polynomial{
public:
    static const int A = 10;

private:
    vector<double> basis;

public:
    polynomial(vector<double> initialVal) : basis(initialVal) {}
    ~polynomial() {}
    vector<double> getBasis() { return basis; }
    polynomial operator+(polynomial);
    polynomial operator*(polynomial);
    polynomial operator*(double);
    polynomial differentiate();
    double value(double x);
    void Change(vector<double> Vec) { basis = Vec; }
    void visualize() const;
    double NewtonIter(double);

};

double polynomial::value(double x){
    double Val{0};
    for (vector<double>::size_type s = 0; s < (*this).basis.size(); s++){
        Val += (*this).basis[s] * pow(x, (double)s);
    }
    return Val;
}

polynomial polynomial::differentiate(){
    polynomial Temp = (*this);
    for (vector<double>::size_type s = 0; s < Temp.basis.size(); s ++){
        Temp.basis[s] *= s;
    }
    for (vector<double>::size_type s = 0; s < Temp.basis.size() - 1; s ++){
        Temp.basis[s] = Temp.basis[s + 1];
    }
    Temp.basis.pop_back();
    return Temp;
}

void polynomial::visualize() const{
    for (vector<double>::size_type s = 0; s < (*this).basis.size(); s ++ ){
        printf("%fx^%d  ", (*this).basis[s], (int)s);
    }
    return;
}

double polynomial::NewtonIter(double x){
    polynomial Poly = *this;
    for (int i = 0; i < 30; i++){
        double y = Poly.value(x);
        polynomial diff = Poly.differentiate();
        double k = diff.value(x);
        if (k == 0.0){
            break;
        }
        x = x - y / k;
    }
    return x;
}

polynomial polynomial::operator+(polynomial P){
    auto self_size = (*this).basis.size(), P_size = P.basis.size();
    if (self_size >= P_size){
        for (vector<double>::size_type s = 0; s < P_size; s++){
            (*this).basis[s] += P.basis[s];
        }
        return (*this);
    }
    else{
        for (vector<double>::size_type s = 0; s < self_size; s++){
            P.basis[s] += (*this).basis[s];
        }
        return P;
    }
}

polynomial polynomial::operator*(polynomial P){
    vector<double> TempVec((*this).basis.size() + P.basis.size() - 1, init);
    polynomial Temp(TempVec);
    for (vector<double>::size_type s = 0; s < (*this).basis.size(); s++){
        for (vector<double>::size_type t = 0; t < P.basis.size(); t++){
            Temp.basis[s + t] += (*this).basis[s] * P.basis[t];
        }
    }
    return Temp;
}

polynomial polynomial::operator*(double x){
    polynomial Temp = (*this);
    for (auto &s : Temp.basis){
        s *= x;
    }
    return Temp;
}