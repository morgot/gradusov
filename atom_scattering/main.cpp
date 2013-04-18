#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <fstream>
using namespace std;

int N, l;
double h, V_0, a, x_max;

double F (double x, double k);
double potential_V ( double x );
double J ( double x, int num );
double j ( double x, int num );
void linSolve(
		  vector<complex<double> >  Dl
		, vector<complex<double> >  D
		, vector<complex<double> >  Du
		, vector<complex<double> >  f
		, vector<complex<double> >  &x
		);
complex<double> amplitude ( double k, vector<complex<double> > y);
complex<double> delta( double k, complex<double> A );
double delta_second( double k, complex<double> A );

int main() {

//===============================================
//=================== INPUT =====================
//===============================================
	
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

//===============================================
//=================== SOVLE =====================
//===============================================
	for( double k = 0.05; k <= 20; k+=0.05 ) {
		vector<complex<double> > Dl(N-1), D(N), Du(N-1), f(N), y(N);
		h = x_max/N;
		double x = h;
		
		//===============================
		 D[0] = complex<double>( -2.0 - 10.0*(h*h/12.0)*F(x,k) ,0);
		Du[0] = complex<double>( 1.0 - (h*h/12.0)*F(x+h,k) ,0);
		
		for ( int i = 1; i < N-1; i++ ) {
			x = (i+1)*h;
			
			Dl[i-1] = complex<double>( 1.0 - (h*h/12.0)*F(x-h,k) ,0);
			D[i] = complex<double>( -2.0 - 10.0*(h*h/12.0)*F(x,k) ,0);
			Du[i] = complex<double>( 1.0 - (h*h/12.0)*F(x+h,k) ,0);
		}
		
		x = N*h;
		Dl[N-2] = complex<double>( 1.0 - (h*h/12.0)*F(x-h,k) ,0);
		D[N-1] = complex<double>( -2.0 - 10.0*(h*h/12.0)*F(x,k) + ( 1.0 - (h*h/12.0)*F(x+h,k))*cos(k*h)
					,( 1.0 - (h*h/12.0)*F(x+h,k))*sin(k*h) );
		//===============================
		
		for ( int i = 0; i < N; i++ ) {
 			x = (i+1)*h;
			
			f[i] = (
			 potential_V(x+h) * J( k*(x+h),l ) +
		      10*potential_V( x ) * J( k*( x ),l ) +
			 potential_V(x-h) * J( k*(x-h),l )
			)*h*h/12.0;
		}
		
		//===============================
		
		linSolve(Dl, D, Du, f, y);
		
		//===============================
		
		for ( int i = 0; i < N; i++ )
			y[i] += J(k*(i+1)*h, l);
		
		//===============================
		
		cout << "=====" << endl;
		for ( int i = 0; i < N; i++ )
			cout << f[i] << endl;
		cout << "=====" << endl;
		
		A = amplitude (k,y);
		
		output << delta_second( k, A );		
		output << endl;
		
		//===============================
		
		vector<complex<double> > Ol(2), O(3), Ou(2), g(3), z(3);
		Ol[0] = complex<double>(0,3);
		Ol[1] = complex<double>(0,6);
		
		O[0] = complex<double>(0,1);
		O[1] = complex<double>(0,4);
		O[2] = complex<double>(0,7);
		
		Ou[0] = complex<double>(0,2);
		Ou[1] = complex<double>(0,5);
		
		g[0] = complex<double>(0,5);
		g[1] = complex<double>(0,26);
		g[2] = complex<double>(0,33);
		
		linSolve(Ol, O, Ou, g, z);
		
		for ( int i = 0; i<3; i++ )
			cout << z[i] << endl;
	}
	
	output.close();
}

double F (double x, double k) {
	return (l*(l+1))/(x*x) + potential_V(x) - double(k*k);
}

double potential_V ( double x ) {
	return (x>=a)||(x<=0) ? 0 : -V_0;
}

double J ( double x, int num ) {
	return x*j(x,num);
}

double j ( double x, int num ) {
	if ( num == 0) {
		if ( x == 0 ) return 1;
		else return (sin(x)/x);
	}
	else if ( num == 1) {
		if (x == 0) return 0;
		return (sin(x)/(x*x)-cos(x)/x);
	}
	else {
		if (x == 0 ) return 0;
		return (2*num + 1)*j(x,num-1)/x - j(x,num-2);
	}
}

void linSolve(
		  vector<complex<double> >  Dl
		, vector<complex<double> >  D
		, vector<complex<double> >  Du
		, vector<complex<double> >  f
		, vector<complex<double> >  &x
		) {
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

complex<double> amplitude ( double k, vector<complex<double> > y) {
	complex<double> S=0;
	double x = h;
	int i = 0;
	while ( x < a ) {
		S += (h/3)*( 
				J(k*(x       ),l) * potential_V(x       ) * y[i  ] +
			      4*J(k*(x + h   ),l) * potential_V(x + h   ) * y[i+1] +
				J(k*(x + 2*h ),l) * potential_V(x + 2*h ) * y[i+2] 
			);
		i+=2;
		x+=2*h;
	}
	
	S /= (k*k);
	return S;
}

complex<double> delta( double k, complex<double> A ){
    return std::log(complex<double>(1,0)+complex<double>(2*k,0)*A*complex<double>(0,1) )/complex<double>(0,2);
}

double delta_second ( double k, complex<double> A ) {
	return atan( 2*k*A.real()/(1-2*k*A.imag()) )/2;
}