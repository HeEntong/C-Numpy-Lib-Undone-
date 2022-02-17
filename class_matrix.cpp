#include <vector>
#include <cstdio>
#include <cmath>
#include "class_polynomial.cpp"
#include <iostream>
#include <algorithm>
using namespace std;
#define init 0.0000000
inline double innerProd(vector<double>, vector<double>);
vector<double> colVec(vector<vector<double>>, double);
auto VecToPoly(vector<vector<double>>) -> vector<vector<polynomial>>;

class matrix{
private:
    vector<vector<double>> elements;
    vector<vector<polynomial>> polyElements;

public:
    matrix() = default;
    matrix(vector<vector<double>> ele) : elements(ele), polyElements(VecToPoly(ele)) {}
    
    vector<vector<double>> getElements() { 
    return elements; }

    vector<vector<polynomial>> getPoly(){
        return polyElements;
    }

    void changeElements(double x, int i, int j){
        (*this).elements[i][j] = x;
        return;
    }

    void changePoly(polynomial x, int i, int j){
        (*this).polyElements[i][j] = x;
        return;
    }


    void add(matrix);

    matrix multiplication(matrix);
    void multiplication(double);

    void print();

    void polyPrint();

    void transpose();

    double Eigen_PowIter(matrix);

    matrix cofactor(int, int);

    double determinant();

    matrix inverse();

    polynomial polyDeterminant();

    matrix polyCofactor(int , int);

    vector<double> Eigen_NewtonIter();
};

void matrix::add(matrix mat){
    for (vector<vector<double>>::size_type s = 0; s < (*this).elements.size(); s++){
        for (vector<double>::size_type t = 0; t < (*this).elements[s].size(); t++){
            (*this).elements[s][t] += (mat).elements[s][t];
        }
    }
}

inline double innerProd(vector<double> V1, vector<double> V2){
    double sum = 0;
    if (!(V1.size() - V2.size())){
        for (vector<double>::size_type s = 0; s < V1.size(); s++){
            sum += V1[s] * V2[s];
        }
    }
    return sum;
}

vector<double> colVec(vector<vector<double>> Vec, double index){
    vector<double> COL;
    for (auto s : Vec){
        COL.push_back(s[index]);
    }
    return COL;
}

matrix matrix::multiplication(matrix mat){
    vector<vector<double>> Result;
    for (vector<vector<double>>::size_type s = 0; s < (*this).elements.size(); s++){
        vector<double> Row;
        for (vector<double>::size_type i = 0; i < (mat).elements[0].size(); i++){
            double S = innerProd((*this).elements[s], colVec((mat).elements, i));
            Row.push_back(S);
        }
        Result.push_back(Row);
    }
    matrix Res(Result);
    return Res;
}

void matrix::multiplication(double coef){
    for (auto &s : (*this).elements){
        for (auto &t : s){
            t *= coef;
        }
    }
    return;
}
void matrix::transpose(){
    vector<double> Row((*this).elements.size(), init);
    vector<vector<double>> TrMat((*this).elements[0].size(), Row);
    for (vector<vector<double>>::size_type s = 0; s < (*this).elements[0].size(); s++){
        for (vector<vector<double>>::size_type t = 0; t < (*this).elements.size(); t++){
            TrMat[s][t] = (*this).elements[t][s];
        }
    }
    (*this).elements = TrMat;
    return;
}

double matrix::Eigen_PowIter(matrix unitVec){
    matrix leftMat = (*this);
    matrix Eigen = (leftMat).multiplication(unitVec);
    for (int i = 0; i < 100; i++){
        Eigen = leftMat.multiplication(Eigen);
        Eigen.multiplication(1 / sqrt
        (innerProd(colVec(Eigen.elements, 0), colVec(Eigen.elements, 0))));
    }
    matrix transMat = Eigen;
    transMat.transpose();
    double lambda = 
    (transMat.multiplication(leftMat)).multiplication(Eigen).elements[0][0] /
     (transMat.multiplication(Eigen)).elements[0][0];
    return lambda;
}

matrix matrix::cofactor(int i, int j){
    matrix Copy = (*this);
    auto pt = Copy.elements.begin();
    pt += i;
    Copy.elements.erase(pt);
    for (auto &s : Copy.elements){
        for (vector<double>::size_type p = j + 1; p < s.size(); p++){
            s[p - 1] = s[p];
        }
        s.pop_back();
    }
    return Copy;
}

matrix matrix::polyCofactor(int i, int j){
    matrix Copy = (*this);
    auto pt = Copy.polyElements.begin();
    pt += i;
    Copy.polyElements.erase(pt);
    for (auto &s : Copy.polyElements){
        for (vector<polynomial>::size_type p = j + 1; p < s.size(); p++){
            s[p - 1] = s[p];
        }
        s.pop_back();
    }
    return Copy;
}

double matrix::determinant(){
    if (!(this->elements.size() - 1)){
        return this->elements[0][0];
    }
    else{
        double DET = {init};
        for (vector<vector<double>>::size_type s = 0; s < this->elements.size(); s++){
            DET += ((s % 2 ?  -1 : 1) * this->elements[s][0])
             * (*this).cofactor(s, 0).determinant() ;
        }
        return DET;
    }
}

matrix matrix::inverse(){
    matrix Copy = (*this);
    double DET = (*this).determinant();
    vector<vector<double>> Storage;
    for (vector<vector<double>>::size_type s = 0; s < Copy.elements.size(); s++){
        vector<double> Temp;
        for (vector<double>::size_type t = 0; t < Copy.elements[0].size(); t++){
            Temp.push_back(Copy.cofactor(t, s).determinant() / DET);
        }
        Storage.push_back(Temp);
    }
    Copy.elements = Storage;
    return Copy;
}

void matrix::print(){
    for (auto i : (*this).elements){
        for (auto j : i){
            printf("%f ", j);
        }
        printf("\n");
    }
}

void matrix::polyPrint(){
    for (auto i : (*this).polyElements){
        for (auto j : i){
            j.visualize();
        }
        printf("\n");
    }
}

auto VecToPoly(vector<vector<double>> VEC) -> vector<vector<polynomial>> {
    vector<vector<polynomial>> Temp_2;
    for (auto &s : VEC){
        vector<polynomial> Temp_1;
        for (auto &t : s){
            polynomial P({t});
            Temp_1.push_back(P);
        }
        Temp_2.push_back(Temp_1);
    }
    return Temp_2;
}

polynomial matrix::polyDeterminant(){
    if ((*this).polyElements[0].size() == 1){
        return polyElements[0][0];
    }
    else{
        /*int Mag = (*this).polyElements[0][0].getBasis().size();
        
        vector<double> Init(Mag, init);*/
        polynomial TempPoly({0.00000});
        for (int  s = 0; s < (int)(*this).getPoly().size(); s++){
            matrix TempCof = (*this).polyCofactor(s, 0);
            polynomial alCof = 
            TempCof.polyDeterminant().multiplication(this->polyElements[s][0]);
            alCof = alCof.multiply(s % 2 ? -1 : 1);
            TempPoly = TempPoly.add(alCof);
        }
        return TempPoly;
        }
}

vector<double> matrix::Eigen_NewtonIter(){
    matrix Temp = *this;
    for (vector<vector<polynomial>>::size_type s = 0; s < Temp.polyElements.size(); s++){
        polynomial K({Temp.polyElements[s][s].getBasis()[0], -1});
        Temp.changePoly(K, s, s);
    }
    polynomial EigenPoly = Temp.polyDeterminant();
    vector<double> EigenVal;
    for (double i = -100.000; i < 100.0000; i++){
        EigenVal.push_back(EigenPoly.NewtonIter(i));
    }
    /*for (auto &s : EigenVal){
        s = ((int)(s * 1000)) / 1000;
    }*/
    sort(EigenVal.begin(), EigenVal.end());
    
    return EigenVal;
}