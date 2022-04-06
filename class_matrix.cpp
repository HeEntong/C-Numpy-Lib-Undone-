#include <vector>
#include <cstdio>
#include <cmath>
#include "class_polynomial.cpp"
#include <iostream>
#include <algorithm>

inline double innerProd(vector<double>, vector<double>);
vector<double> colVec(vector<vector<double>>, double);
auto VecToPoly(vector<vector<double>>) -> vector<vector<polynomial>>;
inline bool isAllZero(vector<double>);

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
        elements[i][j] = x;
        return;
    }

    void changePoly(polynomial x, int i, int j){
        polyElements[i][j] = x;
        return;
    }

    void gaussianElimination();

    matrix operator+(matrix);

    matrix operator*(matrix);
    matrix operator *(double);

    void print() const;

    void polyPrint();

    matrix transpose();

    double Eigen_PowIter(matrix);

    matrix cofactor(int, int);

    double determinant();

    matrix inverse();

    polynomial polyDeterminant();

    matrix polyCofactor(int , int);

    vector<double> Eigen_NewtonIter();

    
};

matrix unitmatrix(double, int);

matrix matrix::operator+(matrix mat){
    for (vector<vector<double>>::size_type s = 0; s < elements.size(); s++){
        for (vector<double>::size_type t = 0; t < elements[s].size(); t++){
            elements[s][t] += (mat).elements[s][t];
        }
    }
    return *this;
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

inline bool isAllZero(vector<double> Vec){
    for (auto s : Vec){
        if (s != 0.00){
            return false;
        }
    }
    return true;
}

vector<double> colVec(vector<vector<double>> Vec, double index){
    vector<double> COL;
    for (auto s : Vec){
        COL.push_back(s[index]);
    }
    return COL;
}

matrix matrix::operator*(matrix mat){
    if (elements[0].size() != mat.elements.size()){
        throw;
    }
    vector<vector<double>> Result;
    for (vector<vector<double>>::size_type s = 0; s < elements.size(); s++){
        vector<double> Row;
        for (vector<double>::size_type i = 0; i < (mat).elements[0].size(); i++){
            double S = innerProd(elements[s], colVec((mat).elements, i));
            Row.push_back(S);
        }
        Result.push_back(Row);
    }
    matrix Res(Result);
    return Res;
}

matrix matrix::operator*(double coef){
    for (auto &s : elements){
        for (auto &t : s){
            t *= coef;
        }
    }
    return (*this);
}

void matrix::gaussianElimination(){
    using std::swap;
    size_t xlen = elements.size();
    size_t ylen = elements[0].size();
    matrix resultMatrix = *this;
    int colNum = 0;
    for (size_t k = 0; k < xlen; k++){
        if (elements[k][colNum] == 0.00){
            for (size_t i = k + 1; i < xlen + 1; i++){
                if (i == xlen){
                    colNum++;
                    break;
                }
                if (elements[i][colNum] != 0.00){
                    swap(elements[i], elements[k]);
                    break;
                }
            }
        }
        if (isAllZero(elements[k])){
            break;
        }
        for (size_t i = k + 1; i < xlen; i++){
            double coef = elements[i][colNum] / elements[k][colNum];
            for (size_t j = 0; j < ylen; j++){
                elements[i][j] -= elements[k][j] * coef;
            }
        }
        for (size_t i = 0; i < k; i++){
            double coef = elements[i][colNum] / elements[k][colNum];
            for (size_t j = 0; j < ylen; j++){
                elements[i][j] -= elements[k][j] * coef;
            }
        }
        colNum++;
    }
}

matrix matrix::transpose(){
    vector<double> Row(elements.size(), init);
    vector<vector<double>> TrMat(elements[0].size(), Row);
    for (vector<vector<double>>::size_type s = 0; s < elements[0].size(); s++){
        for (vector<vector<double>>::size_type t = 0; t < elements.size(); t++){
            TrMat[s][t] = elements[t][s];
        }
    }
    return TrMat;
}

double matrix::Eigen_PowIter(matrix unitVec){
    matrix leftMat = (*this);
    matrix Eigen = leftMat * unitVec;
    for (int i = 0; i < 100; i++){
        Eigen = leftMat * Eigen;
        Eigen = Eigen * (1 / sqrt
        (innerProd(colVec(Eigen.elements, 0), colVec(Eigen.elements, 0))));
    }
    matrix transMat = Eigen;
    transMat.transpose();
    double lambda = 
    ((transMat*leftMat)*Eigen).elements[0][0] /
     (transMat*Eigen).elements[0][0];
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
            Temp.push_back((Copy.cofactor(t, s).determinant() / DET) * ((s + t) % 2 ? -1 : 1));
        }
        Storage.push_back(Temp);
    }
    Copy.elements = Storage;
    return Copy;
}

void matrix::print() const{
    printf("[\n");
    for (auto i : elements){
        printf("  [");
        for (auto j : i){
            printf("%f, ", j);
        }
        printf("]\n");
    }
    printf("]");
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
        polynomial TempPoly({0.00000});
        for (int  s = 0; s < (int)(*this).getPoly().size(); s++){
            matrix TempCof = (*this).polyCofactor(s, 0);
            polynomial alCof = 
            TempCof.polyDeterminant() * this->polyElements[s][0];
            alCof = alCof * (s % 2 ? -1 : 1);
            TempPoly = TempPoly + alCof;
        }
        return TempPoly;
        }
}

matrix unitMatrix(double x, int size){
    vector<double> zeros;
    zeros.assign(size, 0.000);
    vector<vector<double>> unitMat;
    for (int i = 0; i < size; i++){
        unitMat.push_back(zeros);
    }
    for (int i = 0; i < size; i++){
        unitMat[i][i] = x;
    }
    return unitMat;
}

vector<double> matrix::Eigen_NewtonIter(){  
    if (elements[0].size() != elements.size()){
        throw;
    }
    matrix Temp = *this;
    for (vector<vector<polynomial>>::size_type s = 0; s < Temp.polyElements.size(); s++){
        std::vector<double> temp = {Temp.polyElements[s][s].getBasis()[0], -1};
        polynomial K(temp);
        Temp.changePoly(K, s, s);
    }
    polynomial EigenPoly = Temp.polyDeterminant();
    vector<double> EigenVal;
    for (double i = -50.000; i < 50.0000; i++){
        EigenVal.push_back(EigenPoly.NewtonIter(i));
    }
    fastSort(EigenVal, 0, EigenVal.size());
    vector<double> returnVal;
    returnVal.push_back(EigenVal[0]);
    for (size_t i = 0; i < EigenVal.size() - 1; i++){
        if (EigenVal[i+1] - EigenVal[i] > 0.001 ){
            returnVal.push_back(EigenVal[i + 1]);
        }
    }
    return returnVal;
}