#include <iostream>

using namespace std;

#include "alglib/ap.h"
#include "alglib/linalg.h"
#include <fstream>
#include <vector>
#include <cmath>

double potential_V( double r, double a, double V_0 );

void print( alglib::real_2d_array A );

double F( double r, int l, double a, double V_0 );

int main()
{
    // Ввод входных параметров:
    int N, l;
    double h, V_0, a, r_max;
    cout << "Введите количество шагов N: ";
    cin >> N;
    cout << endl;

    cout << "Введите V_0: ";
    cin >> V_0;
    cout << endl;

    cout << "Введите a: ";
    cin >> a;
    cout << endl;

    r_max = 2*a;

    cout << "Введите l: ";
    cin >> l;
    cout << endl;

    // Расчет матриц A и B:

    alglib::real_2d_array A, B;
    alglib::real_1d_array E;
    alglib::real_2d_array f;
    E.setlength( N );
    f.setlength( N, N );
    A.setlength( N, N );
    B.setlength( N, N );

    for ( int i; i < A.rows(); i++)
        for ( int j = 0; j < A.cols(); j++ ) A[i][j] = 0;

    h = r_max/double(N);
    double r = h;
    for ( int i = 0; i < A.rows(); i++){

        if( i != 0 ) {
            A[i][i-1] = F( r-h, l, a, V_0 )*h*h/12 - 1;
            B[i][i-1] = h*h/12;
        }
        A[i][i] = 2 + 5*F( r, l, a, V_0 )*h*h/6;
        B[i][i] = 10*h*h/12;

        if( i != A.cols()-1 ){
            A[i][i+1] = F( r + h, l, a, V_0 )*h*h/12 -1;
            B[i][i+1] = h*h/12;
        }
        r += h;
    }


    // Поиск энергий и состояний:

    bool status = alglib::smatrixgevd( A, N, false, B, false, 1, 1, E, f);

    // Выборка по заданному критерию
    vector< vector< double > > solve;
    solve.clear();
    double eps = 10;
    cout << "Энергии связных состояний: "<< endl;

    for (int i = 0; i < E.length(); i++){

        if ( E[i] > 0 ) break;

        int doub_a = int(1.5*a/h);

        if ( abs( f[ int( doub_a ) ][i] - exp( -sqrt( abs( E[i] ) ) * a ) ) < eps ){
            vector< double > temp( f.cols(), 0);
            cout << "E"<< i <<" = " << E[i] << endl;
            for ( int j = 0; j < f.cols(); j++ ){
                temp[j] = f[j][i];
            }
            solve.push_back(temp);
        }
    }
    // Вывод в файл
    ofstream output("output.txt");

    for(int i=0; i < solve.size(); i++){
        output << "#E = "<< E[i] << endl;
        for(int j = 0; j < solve[0].size(); j++){
            output << h*(j+1) << " "<< solve[i][j] << endl;
        };
    output << endl;
    output << endl;
    };

    output.close();


    return 0;
}

// Определение потенциала V(r):

double potential_V ( double r, double a, double V_0 ){

    if( r > a ) return 0;

    else {
        return -V_0;
    }
}

// Определение пКотенциала F(r, l):

double F( double r, int l, double a, double V_0 ){

    return potential_V(r, a, V_0) + ( double(l)*( double( l + 1 ) ) )/r/r;
}

void print( alglib::real_2d_array A ){
    for(int i=0; i < A.rows(); i++){
        for(int j = 0; j < A.cols(); j++){
            std::cout << A[i][j] << ' ';
        };
    std::cout << std::endl;
    };
};
