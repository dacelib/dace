#include <dace/dace.h>
#include <cmath>
#include <fstream>
using namespace std; using namespace DACE;

const int order = 10;

template<typename T> T Nf( T x, T p )
{
    return x - (x*x-p)/(2*x);
}

// Exercise 4.1.3: Naive Newton
void ex4_1_3( )
{
    double tol = 1e-14;             // tolerance
    double p0 = 4; double x0 = 1;   // expansion point and initial guess
    DA p = p0 + DA(1);              // DA parameter
    DA x = x0;                      // DA initial guess
    DA xp;
    int i = 0;

    do
    {
        xp = x;
        x = Nf( xp, p );
        i++;
    } while( (abs(xp-x) > tol) && (i<1000) );
    cout << "Exercise 4.1.3: Naive Newton" << endl << x << endl << sqrt(p)-x
         << "Number of iterations: " << i << endl << endl;
}

// Exercise 4.1.4: complicated parameters
void ex4_1_4( )
{
    double p0 = 0; double x0 = 1;       // x0 must now satisfy f(x0,cos(p0))=0
    DA p = cos( p0+DA(1) );
    DA x = x0;

    int i = 1;
    while( i <= order )
    {
        x = Nf( x, p );
        i *= 2;
    }
    cout << "Exercise 4.1.3: Naive Newton" << endl << x << endl << sqrt(p)-x;
}

// Exercise 4.2.1: Full DA Newton solver
void ex4_2_1( double p0 )
{
    const double tol = 1e-14;
    double x0 = p0/2.0, xp;   // x0 is just some initial guess
    int i = 0;

    // double precision computation => fast
    do{
        xp = x0;
        x0 = Nf( xp, p0 );
        i++;
    } while( (abs(xp-x0) > tol) && (i<1000) );

    // DA computation => slow
    DA p = p0 + DA(1);
    DA x = x0;
    i = 1;
    while( i <= order )
    {
        x = Nf( x, p );
        i *= 2;
    }

    cout << "Exercise 4.2.1: Full DA Newton" << endl << x << endl << sqrt(p)-x;
}

// Exercise 4.2.2 & 4.2.3: Kepler's equation solver
template<typename T> T NKep( T E, T M, T ecc )
{
    return E - (E-ecc*sin(E)-M)/(1-ecc*cos(E));
}

// double precision Kepler solver
double Kepler( double M, double ecc )
{
    const double tol = 1e-14;
    double E0 = M, Ep;
    int i = 0;
    
    do{
        Ep = E0;
        E0 = NKep( Ep, M, ecc );
        i++;
    } while( (abs(Ep-E0) > tol) && (i<1000) );

    return E0;
}

void ex4_2_2( double M0, double ecc0 )
{
    const double tol = 1e-14;

    DA M = M0 + DA(1);
    DA E = Kepler( M0, ecc0 );      // reference solution
    DA ecc = ecc0;                  // keep eccentricity constant (4.2.2)
    //DA ecc = ecc0 + 0.1*DA(2);    // also expand w.r.t. eccentricity (4.2.3)

    int i = 1;
    while( i <= order )
    {
        E = NKep( E, M, ecc );
        i *= 2;
    }

    cout << "Exercise 4.2.2: Expansion of the Anomaly" << endl
         << "Resulting expansion:" << endl << E << endl
         << "Residual error:" << endl << (E-ecc*sin(E)-M) << endl;

    // sample the resulting polynomial over M0+-3 rad
    ofstream f( "ex4_2_2.dat" );
    for( i = -300; i < 300; i++ )
    {
        f << M0+i/100.0 << "    " << E.evalScalar( i/100.0 )
                        << "    " << Kepler( M0+i/100.0, ecc0 ) << endl;
    }
    f.close( );
    // gnuplot command: plot 'ex4_2_2.dat' u ($1*180/pi):($2*180/pi) w l t 'DA', 'ex4_2_2.dat' u ($1*180/pi):($3*180/pi) w l t 'pointwise'
}

int main( void )
{
    DA::init( order, 2 );      // init with maximum computation order

    ex4_1_3( );
    ex4_1_4( );

    ex4_2_1( 9.0 );
    ex4_2_2( 0.0, 0.5 );
}