#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <fstream>
using namespace std;

double potential_V ( double r, double a, double V_0 );
double F( double x, int l, double a, double V_0, double k );

double j( double x, int l);
double J( double x, int l);
void linSolve( vector<complex<double> >  Dl, vector<complex<double> >  D, vector<complex<double> >  Du, vector<complex<double> >  f, vector<complex<double> >  &x);
void print( vector<complex<double> > A );

complex<double> amplitude( double k, int l, double h, double a, double V_0, vector<complex<double> > y );
complex<double> delta( double k, int l, complex<double> A );

int main()
{
    // Ввод входных параметров:
    int N, l, k;
    double h, V_0, a, x_max;

    complex<double> A;
    cout << "Введите количество шагов N: ";
    cin >> N;
    cout << endl;

    cout << "Введите V_0: ";
    cin >> V_0;
    cout << endl;

    cout << "Введите a: ";
    cin >> a;
    cout << endl;

    x_max = 1.5*a;

    cout << "Введите l: ";
    cin >> l;
    cout << endl;

    ofstream output("output.txt");
    for (double k = 0.05; k<5; k+=0.05){

    // Расчет матрицы A:
    vector<complex<double> > Dl(N-1), D(N), Du(N-1), f(N), y(N);
    h = x_max/N;
    double x = h;

    for ( int i = 0; i <= N-1; i++) {

        if( i != 0 ) {
            Dl[i-1] = complex<double>( 1 - (h*h/12.0)*F( x-h, l, a, V_0, k ), 0 );
        }

        D[i] = complex<double>( (-2)-10*F( x, l, a, V_0, k )*h*h/12.0, 0);

        if( i != N-1 ){
            Du[i] = complex<double>(1 - (h*h/12.0)*F( x+h, l, a, V_0, k ), 0);
        }
        x += h;

    }

    D[N-1] =complex<double>( (-2)-10*F( x-h, l, a, V_0, k )*h*h/12.0+cos(k*h)*(1-F(x,l, a, V_0, k)*h*h/12.0), sin(k*h)*(1-F(x,l, a, V_0, k)*h*h/12.0) );

    x = h;
    for( int i = 0; i < N; i++){
        f[i] = (potential_V(x+h, a, V_0)*J(k*(x+h),l)+potential_V(x, a, V_0)*J(k*x,l)+potential_V(x-h, a, V_0)*J(k*(x-h),l))*(-h)*h/12.0;
        x += h;
    }


    linSolve(Dl, D, Du, f, y);
    x=0;
    for( int i = 0; i < N; i++){
        y[i] += J(x, l);
        x+=h;
    }
    A = amplitude(k, l, h, a, V_0, y);


    output << delta( k, l, A );
        output << endl;

}

    output.close();

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
double F( double x, int l, double a, double V_0, double k ){
    return l*(l+1)/x/x + potential_V(x, a, V_0) - double(k*k);
}

void linSolve( vector<complex<double> >  Dl, vector<complex<double> >  D, vector<complex<double> >  Du, vector<complex<double> >  f, vector<complex<double> >  &x){
    x = f;
    for (int i = 0; i < D.size()-1; i++){
        D[i+1]-=Du[i]*Dl[i]/D[i];
        x[i+1] -= x[i]*Dl[i]/D[i];
    }
    for (int i = D.size()-1; i > 0; i--){

        x[i-1] -= x[i]*Du[i-1]/D[i];
    }
    for (int i=0; i < D.size(); i++){
        x[i]/=D[i];
    }
}
void print( vector<complex<double> > A ){
    for(int j = 0; j < A.size(); j++){
            std::cout << A[j] << std::endl;
        };
};


complex<double> amplitude( double k, int l, double h, double a, double V_0, vector<complex<double> > y ){
    complex<double> S=0;
    double x = h;
    int i=0;
    while( x <= a){
        S+=h/3*( J(x-h,l)*potential_V(x-h, a, V_0)*y[i] + 4*J(x,l)*potential_V(x, a, V_0)*y[i+1] + J(x+h,l)*potential_V(x+h, a, V_0)*y[i+2]);
        x += h;
        i++;
    }
    S /= (k*k);
    return S;
}
complex<double> delta( double k, int l, complex<double> A ){
    return std::log(complex<double>(1,0)+complex<double>(2*k,0)*A*complex<double>(0,1) )/complex<double>(0,2);
}
