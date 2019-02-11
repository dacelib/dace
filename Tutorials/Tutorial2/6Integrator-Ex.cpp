#include <dace/dace.h>
#include <cmath>
#include <fstream>
using namespace std; using namespace DACE;

// Exercise 6.1.1: The Mid-point rule integrator
template<typename T> T midpoint( T x0, double t0, double t1, T (*f)(T,double) )
{
	const double hmax = 0.005;
	int steps = ceil( (t1-t0)/hmax );
	double h = (t1-t0)/steps;
    double t = t0;

    T xmid;
	for( int i = 0; i < steps; i++ )
	{
        xmid = x0 + 0.5*h*f( x0, t );
        x0 = x0 + h*f( xmid, t+0.5*h );
		t += h;
	}

    return x0;
}


// Exercise 6.1.2: Model of the (uncontrolled) pendulum
// x = ( theta, theta_dot )
// since the motion for x decouples we ignore it here
template<class T> AlgebraicVector<T> pendulumRHS( AlgebraicVector<T> x, double t )
{
    // pendulum constants
    static const double l = 0.61;       // m
    static const double m = 0.21;       // kg
    static const double M = 0.4926;     // kg
    static const double g = 9.81;       // kg*m/s^2

    AlgebraicVector<T> res(2);

    res[0] = x[1];
    res[1] = ((M+m)*g*sin(x[0])-m*l*sqr(x[1])*sin(x[0])*cos(x[0]))/((M+m)*l+m*l*sqr(cos(x[0])));

    return res;
}

void ex6_1_2( )
{
    AlgebraicVector<DA> xDA(2);
    AlgebraicVector<double> xdb(2);
    double t;
    const double dt = 0.05;          // take a snap shot every 0.1s

    ofstream file( "ex6_1_2.dat" );
    file.precision( 16 );

    xdb[0] = 0.0; xdb[1] = 0.2;     // initial condition (double)
    t = 0;
    for( int i = 0; i < 100; i++ )
    {
        xdb = midpoint( xdb, t, t+dt, pendulumRHS );    // propagate forward for dt seconds
        file << t << "   " << xdb[0] << "   " << xdb[1] << endl;
        t += dt;
    }
    file << endl << endl;

    xDA[0] = 0.0; xDA[1] = 0.2+0.04*DA(1);  // initial condition (DA)
    t = 0;
    for( int i = 0; i < 100; i++ )
    {
        xDA = midpoint( xDA, t, t+dt, pendulumRHS );    // propagate forward for dt seconds
        file << t << "   " << xDA[0].evalScalar(-1.0) << "   " << xDA[0].evalScalar(+1.0);
        file << "   " << xDA[1].evalScalar(-1.0) << "   " << xDA[1].evalScalar(+1.0) << endl;
        t += dt;
    }
    file << endl;
    file.close( );

    cout << "Exercise 6.1.2: Model of the (uncontrolled) pendulum" << endl << xDA << endl;
}


// Exercise 6.1.3: Set propagation
// the right hand side
template<class T> AlgebraicVector<T> f( AlgebraicVector<T> x, double t )
{
    const double alpha = 0.1;
    AlgebraicVector<T> res(2);

    res[0] = -x[1];
    res[1] =  x[0];
    return (1.0 + alpha*x.vnorm())*res;
}

// convenience routine to evaluate and plot
void plot( AlgebraicVector<DA> x, double t, int N, ofstream &file )
{
    AlgebraicVector<double> arg(2), res;
    for( int i = -N; i <= N; i++ )
    {
        arg[0] = (double)i/N;
        for( int j = -N; j <= N; j++ )
        {
            arg[1] = (double)j/N;
            res = x.eval(arg);
            file << t << "    " << res[0] << "    " << res[1] << endl;
        }
        file << endl;
    }
    file << endl;
}

void ex6_1_3( )
{
    const double pi =  3.141592653589793;
    AlgebraicVector<DA> x(2);
    double t;
    const double dt = 2.0*pi/6.0;

    ofstream file( "ex6_1_3.dat" );
    file.precision( 16 );

    x[0] = 2.0 + DA(1); x[1] = DA(2);  // initial condition box
    t = 0;
    plot( x, t, 7, file );

    for( int i = 0; i < 6; i++ )
    {
        x = midpoint( x, t, t+dt, f );    // propagate forward for dt seconds
        t += dt;
        plot( x, t, 7, file );
    }
    file.close( );

    cout << "Exercise 6.1.3: Set propagation" << endl << x << endl;
}


// Exercise 6.1.4: State Transition Matrix
void ex6_1_4( )
{
    const double pi =  3.141592653589793;
    AlgebraicVector<DA> x(2);

    x[0] = 1.0 + DA(1); x[1] = 1.0 + DA(2);  // initial condition around (1,1)
    x = midpoint( x, 0, 2*pi, f );

    cout << "Exercise 6.1.4: State Transition Matrix" << endl;
    cout.precision( 7 );
    cout << x[0].deriv(1).cons() << "    " << x[0].deriv(2).cons() << endl;
    cout << x[1].deriv(1).cons() << "    " << x[1].deriv(2).cons() << endl;
    cout << endl;
}


// Exercise 6.1.5: Parameter dependence
// the right hand side (note: now it can only be evaluated with DA because alpha is a DA!)
AlgebraicVector<DA> fParam( AlgebraicVector<DA> x, double t )
{
    const DA alpha = 0.05 + 0.05*DA(1);     // parameter, now it's a DA
    AlgebraicVector<DA> res(2);

    res[0] = -x[1];
    res[1] =  x[0];
    return (1.0 + alpha*x.vnorm())*res;
}

void ex6_1_5( )
{
    const double pi =  3.141592653589793;
    AlgebraicVector<DA> x(2);

    x[0] = 1.0; x[1] = 1.0;         // initial condition (1,1)
    x = midpoint( x, 0, 2*pi, fParam );

    ofstream file( "ex6_1_5.dat" );
    file.precision( 16 );
    file << "1 1" << endl << endl << endl;
    for( int i = 0; i <= 20; i++ )
        file << x[0].evalScalar(-1.0+i/10.0) << "   " << x[1].evalScalar(-1.0+i/10.0) << endl;
    file.close( );

    cout << "Exercise 6.1.5: Parameter dependence" << endl << x;
    cout << endl;
}


// Exercise 6.2.1: 3/8 rule RK4 integrator
template<typename T> T rk4( T x0, double t0, double t1, T (*f)(T,double) )
{
	const double hmax = 0.01;
	int steps = ceil( (t1-t0)/hmax );
	double h = (t1-t0)/steps;
    double t = t0;

    T k1, k2, k3, k4;
	for( int i = 0; i < steps; i++ )
	{
        k1 = f( x0, t );
        k2 = f( x0 + h*k1/3.0, t + h/3.0 );
        k3 = f( x0 + h*(-k1/3.0 + k2), t + 2.0*h/3.0 );
        k4 = f( x0 + h*(k1 - k2 + k3), t + h );
        x0 = x0 + h*(k1 + 3*k2 + 3*k3 +k4)/8.0;
		t += h;
	}

    return x0;
}


// Exercise 6.2.2: Artsy Set Propagation
void ex6_2_2( )
{
    const double pi =  3.141592653589793;
    AlgebraicVector<DA> x(2);
    double t;
    const double dt = 2.0*pi/6.0;

    ofstream file( "ex6_2_2.dat" );
    file.precision( 16 );

    // initial condition (c.f. 5Vectors-Ex.cpp)
    x[0] = DA(1);
    x[1] = ((1.0-DA(1)*DA(1))*(DA(2)+1.0)+(DA(1)*DA(1)*DA(1)-DA(1))*(1.0-DA(2)))/2.0;
    t = 0;
    plot( x, t, 7, file );

    for( int i = 0; i < 6; i++ )
    {
        x = midpoint( x, t, t+dt, f );    // propagate forward for dt seconds
        t += dt;
        plot( x, t, 7, file );
    }
    file.close( );

    cout << "Exercise 6.2.2: Artsy Set propagation" << endl << x << endl;
}


// Exercise 6.2.3: CR3BP
template<typename T> AlgebraicVector<T> CR3BP( AlgebraicVector<T> x, double t )
{
    const double MU = 0.30404234e-5;
    AlgebraicVector<T> res(6);
    T d1, d2;

    d1 = sqrt( sqr(x[0]+MU) + sqr(x[1]) + sqr(x[2]) );
    d1 = 1.0/(d1*d1*d1);    // first distance
    d2 = sqrt( sqr(x[0]+MU-1.0) + sqr(x[1]) + sqr(x[2]) );
    d2 = 1.0/(d2*d2*d2);    // second distance

    res[0] = x[3];
    res[1] = x[4];
    res[2] = x[5];
    res[3] = x[0] + 2.0*x[4] - d1*(1-MU)*(x[0]+MU) - d2*MU*(x[0]+MU-1.0);
    res[4] = x[1] - 2.0*x[3] - d1*(1-MU)*x[1]      - d2*MU*x[1];
    res[5] =                 - d1*(1-MU)*x[2]      - d2*MU*x[2];

    return res;
}

void ex6_2_3( )
{
    double T = 3.05923;
    AlgebraicVector<double> x0(6);
    x0[0] = 0.9888426847; x0[1] = 0; x0[2] = 0.0011210277; x0[3] = 0; x0[4] = 0.0090335498; x0[5] = 0;
    AlgebraicVector<DA> x = x0 + AlgebraicVector<DA>::identity( );

    DA::pushTO( 1 );    // only first order computation needed

    x = rk4( x, 0, T, CR3BP );

    cout << "Exercise 6.2.3: CR3BP STM" << endl;
    cout.precision( 6 );
    cout << cons(x);
    for( int i = 0; i < 6; i++ )
    {
        for( int j = 1; j <= 6; j++ )
        {
            cout << cons(x[i].deriv(j)) << "  ";
        }
        cout << endl;
    }
    cout << endl;

    DA::popTO( );
}


// Exercise 6.2.4: Set propagation revisited
void ex6_2_4( )
{
    const double pi =  3.141592653589793;
    AlgebraicVector<DA> x(2);
    double t;
    const double dt = 2.0*pi/6.0;

    ofstream file( "ex6_2_4.dat" );
    file.precision( 16 );

    // initial condition box, in polar coordinates
    x[0] = cos(0.3*DA(2))*(2.0 + DA(1)); x[1] = sin(0.3*DA(2))*(2.0 + DA(1));
    t = 0;
    plot( x, t, 40, file );

    for( int i = 0; i < 6; i++ )
    {
        x = midpoint( x, t, t+dt, f );    // propagate forward for dt seconds
        t += dt;
        plot( x, t, 40, file );
    }
    file.close( );

    cout << "Exercise 6.2.4: Set propagation revisited" << endl << x << endl;
}


// Exercise 6.2.5: The State Transition Matrix reloaded
void ex6_2_5( )
{
    const double pi =  3.141592653589793;
    AlgebraicVector<DA> x(2);

    x[0] = 1.0 + DA(2); x[1] = 1.0 + DA(3);         // initial condition (1,1) plus DA identity (but in DA(2) and DA(3) as DA(1) is already used for alpha!)
    x = midpoint( x, 0, 2*pi, fParam );

    AlgebraicVector<DA> arg(3);
    arg[0] = DA(1);         // we want to evaluate the derivatives at (alpha,0,0), so keep DA(1) and replace DA(2) and DA(3) by zero
    arg[1] = 0;
    arg[2] = 0;

    cout << "Exercise 6.2.5: The State Transition Matrix reloaded" << endl;
    cout << x[0].deriv(2).eval(arg) << x[0].deriv(3).eval(arg) << endl;
    cout << x[1].deriv(2).eval(arg) << x[1].deriv(3).eval(arg) << endl;
    cout << endl;
}


int main( void )
{
    DA::init( 15, 6 );      // init with maximum computation order

    ex6_1_2( );
    ex6_1_3( );
    ex6_1_4( );
    ex6_1_5( );

    ex6_2_2( );
    ex6_2_3( );
    ex6_2_4( );
    ex6_2_5( );
}

