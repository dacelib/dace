#include <dace/dace.h>
#include <cmath>
using namespace std; using namespace DACE;

const double pi =  3.141592653589793;

// Exercise 2.1.1: derivatives
template<typename T> T somb( T x, T y )
{
//    return sin(y)*exp(x*x);
    T r = sqrt( x*x + y*y );
    return sin(r)/r;
}

void ex2_1_1( )
{
    double x0 = 2.0, y0 = 3.0;
    DA x = DA(1); DA y = DA(2);
    DA func = somb( x0+x, y0+y );    // expand sombrero function around (x0,y0)

    // compute the derivative using DA
    double dadx = cons(func.deriv(1));
    double dady = cons(func.deriv(2));
    double dadxx = cons(func.deriv(1).deriv(1));
    double dadxy = cons(func.deriv(1).deriv(2));
    double dadyy = cons(func.deriv(2).deriv(2));
    double dadxxx = cons(func.deriv(1).deriv(1).deriv(1));

    // compute the derivatives using divided differences
    const double h = 1e-3;
    double dx = (somb(x0+h,y0)-somb(x0-h,y0))/(2.0*h);
    double dy = (somb(x0,y0+h)-somb(x0,y0-h))/(2.0*h);
    double dxx = (somb(x0+2.0*h,y0)-2.0*somb(x0,y0)+somb(x0-2.0*h,y0))/(4.0*h*h);
    double dxy = (somb(x0+h,y0+h)-somb(x0-h,y0+h)-somb(x0+h,y0-h)+somb(x0-h,y0-h))/(4.0*h*h);
    double dyy = (somb(x0,y0+2.0*h)-2.0*somb(x0,y0)+somb(x0,y0-2.0*h))/(4.0*h*h);
    double dxxx = (somb(x0+3.0*h,y0)-3.0*somb(x0+h,y0)+3.0*somb(x0-h,y0)-somb(x0-3.0*h,y0))/(8.0*h*h*h);

    cout << "Exercise 2.1.1: Numerical derivatives\n"
         << "d/dx:    " << abs(dadx-dx) << endl
         << "d/dy:    " << abs(dady-dy) << endl
         << "d/dxx:   " << abs(dadxx-dxx) << endl
         << "d/dxy:   " << abs(dadxy-dxy) << endl
         << "d/dyy:   " << abs(dadyy-dyy) << endl
         << "d/dxxx:  " << abs(dadxxx-dxxx) << endl << endl;
}

// Exercise 2.1.2: indefinite integral
void ex2_1_2( )
{
    DA x = DA(1);
    DA func = (1.0/(1+sqr(x))).integ(1);    // DA integral
    DA integral = atan(x);                // analytical integral DA expanded
    cout << "Exercise 2.1.2: Indefinite integral\n" << func-integral << endl;
}

// Exercise 2.1.3: expand the error function
void ex2_1_3( )
{
    DA t = DA(1);
    DA erf = 2.0/sqrt(pi)*exp(-sqr(t)).integ(1);    // error function erf(x)
    cout << "Exercise 2.1.3: Error function\n" << erf << endl;
}

// Exercise 2.2.1: DA based Newton solver
template<typename T> T f( T x )
{
    return x*sin(x)+cos(x);
//    const double e = 0.1; const double M = 1.0;
//    return x-e*sin(x)-M;      // Kepler's equation
}

double ex2_2_1( double x0 )
{
    DA::pushTO( 1 );     // for this Newton solver we only need first derivatives

    const double err = 1e-14;
    DA x = DA(1), func;
    int i = 0;

    do {
        func = f( x0+x );    // expand f around x0
        x0 = x0 - cons(func)/cons(func.deriv(1));       // Newton step
        i++;
    } while( (abs(cons(func))>err) && (i<1000) );

    cout << "Exercise 2.2.1: DA Newton solver\n"
         << "root x0:           " << x0 << endl
         << "residue at f(x0):  " << abs(f(x0)) << endl
         << "Newton iterations: " << i << endl << endl;

    DA::popTO( );       // don't forget to reset computation order to old value for following computations
    return x0;
}

// Exercise 2.2.2: expand the error function around x0
void ex2_2_2( double x0 )
{
    DA t = x0+DA(1);
    DA erf = 2.0/sqrt(pi)*exp(-sqr(t)).integ(1);    // error function erf(x)
    cout << "Exercise 2.2.2: Shifted indefinite integral\n" << erf << endl;
}

int main( void )
{
    DA::init( 30, 2 );
    cout.precision( 15 );

    ex2_1_1( );
    ex2_1_2( );
    ex2_1_3( );

    ex2_2_1( 3.6 );
    ex2_2_2( 1.0 );
}