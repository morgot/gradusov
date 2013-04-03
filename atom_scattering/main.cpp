#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
using namespace std;

double potential_V ( double r, double a, double V_0 );
double F( double x, int l, double a, double V_0, int k );

double j( double x, int l);
double J( double x, int l);

int main()
{
    // Ввод входных параметров:
    int N, l, k;
    double h, V_0, a, x_max;
    cout << "Введите количество шагов N: ";
    cin >> N;
    cout << endl;

    cout << "Введите V_0: ";
    cin >> V_0;
    cout << endl;

    cout << "Введите a: ";
    cin >> a;
    cout << endl;

    x_max = 2*a;

    cout << "Введите l: ";
    cin >> l;
    cout << endl;

    cout << "Введите k: ";
    cin >> k;
    cout << endl;
    // Расчет матрицы A:


    vector<complex<double> > Dl(N-1), D(N), Du(N-1), f(N);
    h = x_max/N;
    double x = h;
    for ( int i = 0; i <= N-1; i++){

        if( i != 0 ) {
            Dl[i-1] = complex<double>( 1 - h*h/12*F( x-h, l, a, V_0, k ), 0 );
        }
        D[i] = complex<double>( (-2)-10*F( x, l, a, V_0, k )*h*h/12, 0);

        if( i != N-1 ){
            Du[i] = complex<double>(1 - h*h/12*F( x+h, l, a, V_0, k ), 0);
        }
        x += h;

    }
    D[N-1] =complex<double>( (-2)-10*F( x, l, a, V_0, k )*h*h/12+cos(k*h), sin(k*h) );
    x = h;
    for( int i = 0; i < N; i++){
        f[i] = (potential_V(x+h, a, V_0)*J(k*(x+h),l)+potential_V(x, a, V_0)*J(k*x,l)+potential_V(x-h, a, V_0)*J(k*(x-h),l))*(-h)*h/12;
        x += h;
    }

}


double potential_V( double r, double a, double V_0 ){
    if( r > a ) return 0;
    else {
          return -V_0;
    }
}

double j( double x, int l){
    if (l == 0){

        if (x==0) return 1;
        return sin( x )/x;

    } else if (l == 1){

        if (x==0) return 0;
        return sin( x )/x/x - cos( x )/x;

    } else {

        if (x==0) return 0;
        return (2*l-1)*j(x, l-1)/x - j(x, l-2);
    }

}

double J( double x, int l){
    return x*j(x, l);
}
double F( double x, int l, double a, double V_0, int k ){
    return l*(l+1)/x/x + potential_V(x, a, V_0) - double(k*k);
}
