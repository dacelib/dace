#include <dace/dace.h>
#include <cmath>
using namespace std;
using namespace DACE;


template<typename T> AlgebraicVector<T> TBP( AlgebraicVector<T> x, double t )
{
    
    AlgebraicVector<T> pos(3), res(6);
    
    pos[0] = x[0]; pos[1] = x[1]; pos[2] = x[2];
    
    T r = pos.vnorm();
    
    const double mu = 398600; // km^3/s^2
    
    res[0] = x[3];
    res[1] = x[4];
    res[2] = x[5];
    
    res[3] = -mu*pos[0]/(r*r*r);
    res[4] = -mu*pos[1]/(r*r*r);
    res[5] = -mu*pos[2]/(r*r*r);
    
    return res;
    
}


template<typename T> AlgebraicVector<T> euler( AlgebraicVector<T> x, double t0, double t1 )
{
	const double hmax = 0.1;
	int steps = ceil((t1-t0)/hmax);
	double h = (t1-t0)/steps;
    double t = t0;

	for( int i = 0; i < steps; i++ )
	{
        x = x + h*TBP(x,t);
		t += h;
	}
    
    return x;
}


int main( void )
{
    
    DA::init( 3, 6 );
    
    AlgebraicVector<DA> x0(6), xf(6);
    
    // Set initial conditions
    const double mu = 398600;
    const double ecc = 0.5;
    
    x0[0] = 6678.0 + DA(1); // 300 km altitude
    x0[1] = 0.0    + DA(2);
    x0[2] = 0.0    + DA(3);
    x0[3] = 0.0    + DA(4);
    x0[4] = sqrt(mu/6678.0)*sqrt(1+ecc) + DA(5);
    x0[5] = 0.0 + DA(6);

    // integrate for half the orbital period
    
    const double pi = 4.0*atan(1.0);

    double a = 6678.0/(1-ecc);

    xf = euler( x0, 0.0, pi*sqrt(a*a*a/mu));
    
    // print initial and final conditions
    
    cout << endl << "Initial conditions:" << endl << endl;
    cout << x0 << endl << endl;
    
    cout << endl << "Final conditions:" << endl << endl;
    cout << xf << endl << endl;
    
    getchar();
    
    // Evaluate for a displaced initial condition
    
    AlgebraicVector<double> Deltax0(6);
    
    Deltax0[0] =  1.0; // km
    Deltax0[1] = -1.0;
    Deltax0[2] =  0.0;
    Deltax0[3] =  0.0;
    Deltax0[4] =  0.0;
    Deltax0[5] =  0.0;
    
    cout << "Displaced Initial condition:" << endl << endl;
    cout << x0.cons() + Deltax0 << endl << endl;
    
    cout << "Displaced Final condition:" << endl << endl;
    cout << xf.eval(Deltax0) << endl;

}



