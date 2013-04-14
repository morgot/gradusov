#include <iostream>

using namespace std;
#include "alglib_headers/ap.h"
#include "alglib_headers/linalg.h"
#include <cmath>

int main()
{
    // Ввод входных параметров:

    // Расчет матриц A и B:

    const int N = 4;

    const double a[] = { 13.00773, 1.962079, 0.444529, 0.1219492 };

    const double Pi = 3.14159265;

    alglib::real_2d_array A, B;
    alglib::real_1d_array E;
    alglib::real_2d_array f;
    E.setlength( N );
    f.setlength( N, N );
    A.setlength( N, N );
    B.setlength( N, N );
    for ( int i; i < A.rows(); i++)
        for ( int j = 0; j < A.cols(); j++ ) A[i][j] = 0;


    for ( int i = 0; i < A.rows(); i++){
        for ( int j = 0; j < A.cols(); j++){

            A[i][j] = 3 * pow(Pi, 1.5) * a[i] * a[j] / pow( a[i]+a[j] , 2.5 );
            A[i][j] -= 2 * Pi / ( a[i] + a[j] );
            B[i][j] = pow( Pi / ( a[i] + a[j] ) , 1.5 );
        }
    }


        // Поиск энергий и состояний:

    bool status = alglib::smatrixgevd( A, N, false, B, false, 1, 1, E, f);

    cout << "E = " << E[0] << endl;

    return 0;
}

