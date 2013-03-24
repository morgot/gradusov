#include <iostream>

using namespace std;

#include "alglib/ap.h"
#include "alglib/linalg.h"

double potential_V( double r );

void print( alglib::real_2d_array A );

double F( double r, int l );

int main()
{
    // Ввод входных параметров:
    int N, l;
    double h, r_max, E;
    cout << "Введите количество шагов N: ";
    cin >> N;
    cout << endl;

    cout << "Введите r_max: ";
    cin >> r_max;
    cout << endl;

    cout << "Введите l: ";
    cin >> l;
    cout << endl;

    alglib::real_2d_array A, B;
    A.setlength( N, N );
    B.setlength( N, N );

    for ( int i; i < A.rows(); i++)
        for ( int j = 0; j < A.cols(); j++ ) A[i][j] = 0;

    h = r_max/double(N);
    double r = h;
    for ( int i = 0; i < A.rows(); i++){

        if( i != 0 ) {
            A[i][i-1] = F( r-h, l )*h*h/12 - 1;
            B[i][i-1] = h*h/12;
        }
        A[i][i] = 2 + 5*F( r, l )*h*h/6;
        B[i][i] = 10*h*h/12;

        if( i != A.cols()-1 ){
            A[i][i+1] = F( r + h, l)*h*h/12 -1;
            B[i][i+1] = h*h/12;
        }
        r += h;
    }

    print(A);
    cout << endl;
    print(B);



    return 0;
}

// Определение потенциала V(r):

double potential_V ( double r ){

    double k = 2.0;
    return k/r;
}

// Определение потенциала F(r, l):

double F( double r, int l ){

    return potential_V(r) + ( double(l)*( double( l + 1 ) ) )/r/r;
}

void print( alglib::real_2d_array A ){
    for(int i=0; i < A.rows(); i++){
        for(int j = 0; j < A.cols(); j++){
            std::cout << A[i][j] << ' ';
        };
    std::cout << std::endl;
    };
};
