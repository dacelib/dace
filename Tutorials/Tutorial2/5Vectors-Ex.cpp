#include <dace/dace.h>
#include <cmath>
#include <fstream>
using namespace std; using namespace DACE;

const unsigned int order = 20;

// Exercise 5.1.1: tangents and normals
template<typename T> T f( AlgebraicVector<T> x )
{
    return -1/3*(sqr(x[0])+sqr(x[1])/2)+exp(x[0]/2+x[1]/4);
}

void ex5_1_1( )
{
    AlgebraicVector<DA> surf(3);
    surf[0] = DA(1);
    surf[1] = DA(2);
    surf[2] = f( surf );   // trick: f() only uses the first 2 components of surf, which we already set
    AlgebraicVector<DA> t1 = surf.deriv(1);
    AlgebraicVector<DA> t2 = surf.deriv(2);
    AlgebraicVector<DA> n = t1.cross(t2).normalize();   // normalized surface normal

    cout << "Exercise 5.1.1: tangents and normals" << endl
         << t1 << t2 << n;
}

// Exercise 5.1.2: (Uncontrolled) Equations of motion of the inverted pendulum
template<typename T> AlgebraicVector<T> ode_pendulum( AlgebraicVector<T> x/*, T u*/ )
{
    // constants
    static const double l = 1.0;    // length of pendulum (m)
    static const double m = 0.1;    // weight of balanced object (kg)
    static const double M = 0.4;    // weight of cart (kg)
    static const double g = 9.81;   // gravity acceleration constant on earth (kg*m/s^2)

    // variables
    AlgebraicVector<T> rhs(2);  // right hand side of the ODE (2 dimensional)
    T sint = sin(x[0]);         // sine of theta
    T cost = cos(x[0]);         // cosine of theta

    // Equations of motion
    rhs[0] = x[1];
    rhs[1] = (/*u+*/(M+m)*g*sint-m*l*sqr(x[1])*sint*cost)/((M+m)*l+m*l*sqr(cost));

    return rhs;
}

void ex5_1_2( )
{
    AlgebraicVector<double> x(2);
    x[0] = 1.0; x[1] = 0.0;
    cout << "Exercise 5.1.2: Equations of Motion" << endl << ode_pendulum( x );
}

// Exercise 5.2.1: Solar flux
void ex5_2_1( )
{
    AlgebraicVector<DA> surf(3);
    surf[0] = DA(1);
    surf[1] = DA(2);
    surf[2] = f( surf );   // trick: f() only uses the first 2 components of surf, which we already set
    AlgebraicVector<DA> t1 = surf.deriv(1).normalize(); // normalizing these helps keep the coefficents small and prevents roundoff errors
    AlgebraicVector<DA> t2 = surf.deriv(2).normalize();
    AlgebraicVector<DA> n = t1.cross(t2).normalize();   // normalized surface normal
    AlgebraicVector<double> sun(3);                     // sun direction
    sun[0] = 0.0; sun[1] = 0.0; sun[2] = 1.0;
    DA flux = n.dot(sun);                               // solar flux on the surface

    // Output results
    cout << "Exercise 5.2.1: Solar flux" << endl << flux;
    const int N = 30;
    AlgebraicVector<double> arg(2);
    AlgebraicVector<double> res(3);
    ofstream file("ex5_2_1.dat");
    for( int i = -N; i <= N; i++ )
    {
        arg[0] = (double)i/N;
        for( int j = -N; j <= N; j++ )
        {
            arg[1] = (double)j/N;
            res = surf.eval(arg);
            file << res[0] << "    " << res[1] << "    " << f( res ) << "    " << flux.eval( arg ) << endl;
        }
        file << endl;
    }
    file.close( );
}

// Exercise 5.2.2: Area
void ex5_2_2( )
{
    DA x = DA(1);
    DA t = DA(2);
    AlgebraicVector<DA> res(2);
    res[0] = DA(1);
    res[1] = ((1.0-x*x)*(t+1.0)+(x*x*x-x)*(1.0-t))/2.0;

    cout << "Exercise 5.2.2: Area" << endl << res;
}

int main( void )
{
    DA::init( order, 2 );      // init with maximum computation order

    ex5_1_1( );
    ex5_1_2( );

    ex5_2_1( );
    ex5_2_2( );
}
