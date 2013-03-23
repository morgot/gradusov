/*
 * matrix.cpp
 *
 *  Created on: 19.02.2013
 *      Author: gregory
 */
#include "matrix.h"
//**************************************************************************************
// Размеры матрицы
int Simple_Matrix::colSize() const{
    return M.size();
}

int Simple_Matrix::rawSize() const{
    return M[0].size();
}
// Конструкторы
Simple_Matrix::Simple_Matrix( int m, int n ) : M( m , std::vector<double>(n, 0)){

}

Simple_Matrix::Simple_Matrix(const Simple_Matrix& N) : M( N.colSize(), std::vector<double>(N.rawSize(), 0) ){
    for(int i=0; i < M.size(); i++){
        for(int j = 0; j < M[i].size(); j++){
            M[i][j] = N[i][j];
        };
    };
}

Simple_Matrix::Simple_Matrix(const std::vector<double>& N) : M( N.size(), std::vector<double>(1, 0) ){

    for ( int i = 0; i < M.size(); i++ ){
        for( int j = 0; j < M[i].size(); j++ ){
            M[i][j] = N[i];
        };
    };
}

// Индексация
const std::vector<double>& Simple_Matrix::operator[]( int i ) const {
    return M[ i ];
}

std::vector<double>& Simple_Matrix::operator[]( int i ){
return M[ i ];
}

// Сложение
Simple_Matrix& Simple_Matrix::operator +=(Simple_Matrix& N){

    for (int i=0; i < M.size(); i++){
        for (int  j = 0; j < M[i].size(); j++){
            M[i][j] += N[i][j];
        };
    };
    return *this;
}

Simple_Matrix operator+( Simple_Matrix A, Simple_Matrix B ){
    Simple_Matrix C = A;
    return C += B;
}

std::vector<double> operator+( std::vector<double> A, std::vector<double> B ){
    std::vector<double> C = A;
    for ( int i = 0; i < C.size(); i++ ){
        C[i] += B[i];
    };
    return C;
}
// Вычитание

Simple_Matrix& Simple_Matrix::operator -=(Simple_Matrix& N){

    for (int i=0; i < M.size(); i++){
        for (int  j = 0; j < M[i].size(); j++){
            M[i][j] -= N[i][j];
        };
    };
    return *this;
}

Simple_Matrix operator-( Simple_Matrix A, Simple_Matrix B ){
    Simple_Matrix C = A;
    return C -= B;
}


std::vector<double> operator-( std::vector<double> A, std::vector<double> B ){
    std::vector<double> C = A;
    for ( int i = 0; i < C.size(); i++ ){
        C[i] -= B[i];
    };
    return C;
}

// Умножение
Simple_Matrix& Simple_Matrix::operator *=( double a ){

    for (int i=0; i < M.size(); i++){
        for (int  j = 0; j < M[i].size(); j++){
            M[i][j] *= a;
        };
    };
    return *this;
};

Simple_Matrix& Simple_Matrix::operator*=( Simple_Matrix& N){

    Simple_Matrix C( this->colSize(), N.rawSize() );
    for( int i = 0; i < C.colSize(); i++ ){
        for( int j = 0; j < C.rawSize(); j++ ){
            for( int k=0; k<this->rawSize(); k++ ){
                C[i][j]+=(*this)[i][k]*N[k][j];
            };
        };
    };
    return *this = C;
}

/*
Simple_Matrix& Simple_Matrix::operator*=( std::vector<double>& N){

    Simple_Matrix C( this->colSize(), 1 );
    for( int i = 0; i < C.colSize(); i++ ){
        for( int j = 0; j < C.rawSize(); j++ ){
            for( int k=0; k<this->rawSize(); k++ ){
                C[i][j]+=(*this)[i][k]*N[k];
            };
        };
    };
    return *this = C;
};
*/

Simple_Matrix operator*( Simple_Matrix A, Simple_Matrix B ){
    Simple_Matrix C = A;
    return C *= B;
};




Simple_Matrix operator*( double a, Simple_Matrix B ){
    Simple_Matrix C = B;
    return C *= a;
};

Simple_Matrix operator*( Simple_Matrix A, double b ){
    Simple_Matrix C = A;
    return C *= b;
};

/*
Simple_Matrix operator*( Simple_Matrix A, std::vector<double> b ){
    Simple_Matrix C = A;
    return C *= b;
};
*/

std::vector<double> operator*(std::vector<double> a, double b){
    std::vector<double> c = a;
    for( int i = 0; i < a.size(); i++ ){
        c[i]*=b;
    };
    return c;
};


std::vector<double> operator*( double a, std::vector<double> b){
    std::vector<double> c = b;
    for( int i = 0; i < b.size(); i++ ){
        c[i]*=a;
    };
    return c;
};

// Деление

Simple_Matrix& Simple_Matrix::operator/=( double b){

    for ( int i = 0; i < M.size(); i++ ){
        for (int  j = 0; j < M[i].size(); j++ ){
            M[i][j] /= b;
        };
    };
    return *this;
}

Simple_Matrix operator/( Simple_Matrix A, double b ){
    Simple_Matrix C = A;
    return C /= b;
};


// Вывод

void Simple_Matrix::print(){
    for(int i=0; i < this->colSize() ; i++){
        for(int j = 0; j < this->rawSize(); j++){
            std::cout << M[i][j] << ' ';
        };
    std::cout << std::endl;
    };
};

// Транспонирование

Simple_Matrix& Simple_Matrix::transpose(){

    Simple_Matrix T = *this;
    *this = Simple_Matrix( T.rawSize(), T.colSize() );
    for( int i = 0; i < this->colSize(); i++ ){
        for( int j = 0; j < this->rawSize(); j++ ){

            (*this)[i][j] = T[j][i];
        };
    };

    return *this;
};

Simple_Matrix transpose(Simple_Matrix A){
    Simple_Matrix B = A;
    return B.transpose();
}

// Скалярное произведение


double dotProd( Simple_Matrix A, Simple_Matrix B){
    double sum = 0;

    for (int i = 0; i < A.colSize(); i++ ){
        sum += A[i][0]*B[i][0];
    }
    return sum;
}







//*******************************************************************************
/*

// Размеры матрицы
int Symmetric_Matrix::colSize() const{
    return n;
}

int Symmetric_Matrix::rawSize() const{
    return n;
}
// Конструкторы
Symmetric_Matrix::Symmetric_Matrix( int m ) : M( m*(m+1)/2 , 0 ){
n = m;
}

Symmetric_Matrix::Symmetric_Matrix(const Symmetric_Matrix& N) :  M( (N.colSize()+1)*N.colSize()/2, 0 ){
    if ( N.colSize() == N.rawSize() ){
        for( int j = 0; j < N.rawSize(); j++ ){
            for( int i = 0; i <= j; i++ ){
                M[ i + j*(j-1)/2 ] = N[ i ] [ j ];
            };
        };
    }
    n = N.rawSize();
}
// Индексация
const std::vector<double>& Symmetric_Matrix::operator[]( int i ) const {

}

std::vector<double>& Symmetric_Matrix::operator[]( int i ){

}

// Сложение
Symmetric_Matrix& Symmetric_Matrix::operator +=(Symmetric_Matrix& N){
}

Symmetric_Matrix operator+( Symmetric_Matrix A, Symmetric_Matrix B ){
}

// Умножение
Symmetric_Matrix& Symmetric_Matrix::operator *=( double a ){
}

Symmetric_Matrix& Symmetric_Matrix::operator*=( Symmetric_Matrix& N){
}

Symmetric_Matrix operator*( Symmetric_Matrix A, Symmetric_Matrix B ){
}
*/
