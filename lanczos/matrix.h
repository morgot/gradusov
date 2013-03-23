/*
 * matrix.h
 *
 *  Created on: 19.02.2013
 *      Author: gregory
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <iostream>

class Simple_Matrix {
    std::vector<std::vector<double> > M;
public:

    Simple_Matrix(int m, int n);
    Simple_Matrix(const Simple_Matrix& N);
    Simple_Matrix(const std::vector<double>& N);

    int colSize() const;
    int rawSize() const;


    const std::vector<double>& operator[]( int ) const;
    std::vector<double>& operator[]( int );

    Simple_Matrix& operator+=( Simple_Matrix& );
    Simple_Matrix& operator -=( Simple_Matrix& );

    Simple_Matrix& operator*=( double );
    Simple_Matrix& operator*=( Simple_Matrix& );

    Simple_Matrix& operator/=( double );
//    Simple_Matrix& operator*=( std::vector<double>& );
    void print();
    Simple_Matrix& transpose();
};

Simple_Matrix operator+( Simple_Matrix A, Simple_Matrix B );
Simple_Matrix operator-( Simple_Matrix A, Simple_Matrix B );
std::vector<double> operator+( std::vector<double> A, std::vector<double> B );
std::vector<double> operator-( std::vector<double> A, std::vector<double> B );

Simple_Matrix operator*( Simple_Matrix A, Simple_Matrix B );

Simple_Matrix operator*( double a, Simple_Matrix B );
Simple_Matrix operator*( Simple_Matrix A, double b );
//Simple_Matrix operator*( Simple_Matrix A, std::vector<double> b );

std::vector<double> operator*(std::vector<double> , double);
std::vector<double> operator*( double , std::vector<double>);

Simple_Matrix operator/( Simple_Matrix A, double b );

double dotProd( Simple_Matrix, Simple_Matrix );
Simple_Matrix transpose(Simple_Matrix A);
#endif /* MATRIX_H_ */
