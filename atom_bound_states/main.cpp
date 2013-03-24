#include <iostream>

using namespace std;

#include "alglib/ap.h"
#include "alglib/linalg.h"

double potential_V( double r );

double F( double r, int l );

int main()
{
    // Ввод входных параметров:
    int N, l;
    double h, r_max, E;
    cout << "Введите количество шагов N: ";
    cin >> n;
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
