/*
 * main.cpp
 *
 *  Created on: 19.02.2013
 *      Author: gregory
 */
#include <iostream>
#include <cmath>
#include "matrix.h"
#include "sym_matrix.h"

bool my_signbit (double num) { return num>0 ? 0 : 1; };

double delta( double lambda, sym_matrix E );

void range( sym_matrix E, double& min, double& max );

void localize( sym_matrix E, double max, double min, double e, std::vector<std::vector<double> >& intervals);

double rootFind( double left, double right, double (*func)( double lambda, sym_matrix E ), double e, double f, sym_matrix E );

int main (){
/*
    int m, n;
    std::cout << "Введите m: ";
    std::cin >> m;
    std::cout << "Введите n: ";
    std::cin >> n;

    std::cout << "Введите матрицу A" << std::endl;
    Simple_Matrix A( n, n ), V( m + 1, n ), E(m,m);
    for (int i=0; i < A.colSize(); i++){
        for (int j = i; j < A.rawSize(); j++){
            std::cin >> A[i][j];
            A[j][i] = A[i][j];
        };
    };
*/

    int m;
  //  std::cout << "Введите n: ";
  //  std::cin >> n;
    std::cout << "Введите m: ";
    std::cin >> m;
    sym_matrix A = read_sym_matrix("matrix.txt");
    A.print();
    int n = A.size();
    sym_matrix E(m);
    Simple_Matrix V( m + 1, n );
    V[0][0] = 1;
    for ( int i = 1; i < V.rawSize(); i++ ){
    V[0][i] = 0;
    }


    double Bj = 0, a = 0;

    for ( int j = 0; j < m; j++){

        Simple_Matrix w(n,1);


        if (j == 0) w = A * V[j]; else w = A * V[j] - Bj*V[j-1];
        a = dotProd( w, V[j] );
        w = w - a * V[j];
        if (!( Bj = std::sqrt(dotProd( w, w ) ) ) ) break;

        w /=Bj;
        for ( int i = 0; i < V[j+1].size(); i++){
            V[j+1][i] = w[i][0];
        }

        E.element( j, j ) = a;
        if (j <= m-2){
            E.element( j, j+1 ) = Bj;
            E.element( j+1, j) = Bj;
        }
    };

    V.transpose();
    E.print();
    /*
    double lambda, e;
    e = 0.001;
    lambda = 5;
    while( std::abs(delta( lambda, E )[0] ) > e){
        lambda = lambda - delta( lambda, E )[0] / delta( lambda, E ) [1];
    }

    std::cout << lambda << std::endl;
*/

    double min, max;

    range( E, min, max );

    double e = 0.001;

    std::vector<std::vector<double> > intervals(0,std::vector<double>(2,0));


    localize( E, max, min, e, intervals );

    for (int i = 0 ; i < intervals.size(); i ++){
        std::cout << intervals[i][1] << ' ' << intervals[i][0] << std::endl;
    }
     std::cout << std::endl;
    double eps = 0.01;
    std::cout << "Интервалы:"<< min << ' ' << max << std::endl;

    for ( int i = 0; i < intervals.size(); i ++){
        std::cout << rootFind( intervals[i][1], intervals[i][0], &delta, eps, delta(intervals[i][1], E), E ) <<" " << delta(rootFind( intervals[i][1], intervals[i][0], &delta, eps, delta(intervals[i][1], E), E ), E)  <<std::endl;
    }
    return 0;
}


double delta( double lambda, sym_matrix E ){

    double d0, d1, t, D0, D1, T;
    double del;

    d1 = 1;
    d0 = 0;
    for ( int i = 0; i < E.size(); i++){
        t = d1;
        if ( i == 0 ){
            d1 = ( E.element( i , i ) - lambda ) * d1;
        } else {
        d1 = ( E.element( i, i ) - lambda ) * d1 - E.element( i-1, i ) * E.element( i-1, i ) * d0;
        }

        d0 = t;
    }
    del = d1;
    return del;
};


void range( sym_matrix E, double& min, double& max ){
    double r;

    for (int i = 0; i < E.size(); i++){
        if ( i == 0){
            r = E.element( i, i+1 );
            min = E.element( i, i ) - r;
            max = E.element( i, i ) + r;
        } else if ( i == E.size()-1 ){
            r = E.element( i, i-1 );
        } else {
            r = E.element( i, i+1 ) + E.element( i, i-1 );
        }
        if ( i != 0 && min > E.element( i, i ) - r) min = E.element( i, i ) - r;
        if ( i != 0 && max < E.element( i, i ) + r) max = E.element( i, i ) + r;

    }
}

void localize( sym_matrix E, double max, double min, double e, std::vector<std::vector<double> >& intervals){

    int m = E.size();
    double step = ( max - min )/m;

    std::vector<double> temp(2,0);

    while ( (step > e) && (intervals.size() < m) ){
        intervals.clear();
        temp[0] = min;
        temp[1] = min;
        bool sign_left, sign_right;

        while ( temp[0] <= max ){
            sign_left = std::signbit( delta( temp[1], E ) );
            sign_right = std::signbit( delta( temp[0], E ) );

            if ( sign_left != sign_right )
                 intervals.push_back( temp );

            temp[1] = temp[0];
            temp[0] += step;
        }

        step /=m;



   }
}

double rootFind( double left, double right, double (*func)( double lambda, sym_matrix E ), double e, double f, sym_matrix E ){
 //   int k;
    double x = (right + left)/2.0;
    if ( std::abs( func( x, E ) ) < e){  //std::abs(std::abs( func( x, E ) ) - std::abs(f) ) < e
        return x;
    }
    double A = func( left, E );
    double B = func( x, E );
    double C = func( right, E );

    f = func( x, E );

    //std::cout <<  func( left, E )  << "    ";
//    std::cout <<  func( x, E ) << "    ";
 //   std::cout <<  func( right, E )  << std::endl;
 //   std::cin >> k;
    bool i = my_signbit( func( x, E ) );
   // std::cout <<  i  << "    "<<  my_signbit( func( left, E ) ) << std::endl;
    if ( i != my_signbit( func( left, E ) ) ){
        return rootFind( left, x, func , e, f, E );
    } else {
        return rootFind( x, right, func , e, f, E );
    }

}
