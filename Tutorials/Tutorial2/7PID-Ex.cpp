#include <dace/dace.h>
#include <cmath>
#include <algorithm>
#include <fstream>
using namespace DACE; using namespace std;

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


// Exercise 7.1.1: Model of the controlled pendulum
// x = ( theta, theta_dot, u )
template<class T> AlgebraicVector<T> pendulumRHS( AlgebraicVector<T> x, double t )
{
    // pendulum constants
    static const double l = 0.61;       // m
    static const double m = 0.21;       // kg
    static const double M = 0.4926;     // kg
    static const double g = 9.81;       // kg*m/s^2

    AlgebraicVector<T> res(3);
    T sint, cost;
    sint = sin(x[0]);
    cost = cos(x[0]);

    res[0] = x[1];
    res[1] = (x[2]+(M+m)*g*sint-m*l*sqr(x[1])*sint*cost)/((M+m)*l+m*l*sqr(cost));
    res[2] = 0;         // u is assumed constant unless changed externally by the controller

    return res;
}


// Exercise 7.1.2: Tuning the PID simulator (double)
void ex7_1_2( )
{
    const double Kp = 8, Ti = 3, Td = 0.3;                  // PID parameters
    const double setPt = 0.0, dt = 0.05;                    // Set point, controller time step (50 ms)
    const double umax = 10;                                 // maximum control (Exercise 7.2.1)

    ofstream file( "ex7_1_2.dat" );

    double intErr = 0, lastErr = 0;                         // PID controller variables

    double t = 0;
    AlgebraicVector<double> x(3);
    x[0] = 0.1; x[1] = 0; x[2] = 0;     // Initial condition (u(0)=x[2]=0)

    // propagate the model for 100 sec
    while( (t<100.0) && (abs(x[0])<1.5) )
    {
        // compute the PID control at this controller time step
        double  err = setPt-x[0];
        double derr = (err-lastErr)/dt;
        intErr += lastErr*dt;
        lastErr = err;
        x[2] = Kp*(err + Td*derr + intErr/Ti);
        // prevent control saturation (Exercise 7.2.1)
//        x[2] = tanh(x[2]/2.0/umax)*umax*2.0;
//        x[2] = max(min(x[2],umax),-umax);

        // output and propagate one time step
        file << t << "   " << x[0] << "   " << x[2] << endl;
        x = rk4( x, t, t+dt, pendulumRHS );
        t += dt;
    }

    cout << "Final angle theta:" << x << endl;
    if(abs(x[0])>1.5)
        cout << "WHOOPSY: Fell over after " << t << " seconds." << endl;

    file.close();
}


// Exercise 7.1.3: PID simulator (DA)
void ex7_1_3( )
{
    const double Kp = 8, Ti = 3, Td = 0.3;                 // PID parameters
    const double setPt = 0.0, dt = 0.05;                    // Set point, controller time step (50 ms)
    const double umax = 10;                                 // maximum control (Exercise 7.2.1)

    ofstream file( "ex7_1_3.dat" );

    DA intErr = 0, lastErr = 0;                             // PID controller variables

    double t = 0;
    AlgebraicVector<DA> x(3);
    x[0] = 0.1+0.1*DA(1); x[1] = 0; x[2] = 0;               // Initial condition

    // propagate the model state for 100 sec
    while( (t<40.0) && (abs(cons(x[0]))<1.5) )
    {
        // compute the PID control
        DA err = setPt-x[0];
        DA derr = (err-lastErr)/dt;
        intErr += lastErr*dt;
        lastErr = err;
        x[2] = Kp*(err + Td*derr + intErr/Ti);
        // prevent control saturation (Exercise 7.2.1)
//        x[2] = tanh(x[2]/umax)*umax;

        // output and propagate one time step (Exercise 7.1.4)
        Interval bx = x[0].bound( );
        Interval bu = x[2].bound( );
        file << t
        << "   " << cons(x[0]) << "   " << bx.m_ub << "   " << bx.m_lb
        << "   " << x[0].evalScalar(-1.0) << "   " << x[0].evalScalar(1.0)
        << "   " << cons(x[2]) << "   " << bu.m_ub << "   " << bu.m_lb
        << "   " << x[2].evalScalar(-1.0) << "   " << x[2].evalScalar(1.0) << endl;
        x = rk4( x, t, t+dt, pendulumRHS );
        t += dt;
    }

    cout << "Final angle theta:" << x << endl;
    if(abs(cons(x[0]))>1.5)
        cout << "WHOOPSY: Fell over after " << t << " seconds." << endl;

    file.close();
}


// Exercise 7.1.5: Bounding
void ex7_1_5( )
{
    DA x = DA(1); DA y = DA(2);
    DA func = sin(x/2)/(2+cos(y/2+x*x));

    Interval a, b, c;

    // bound by rasterizing
    AlgebraicVector<double> arg(2);
    a.m_lb = 9999999; a.m_ub = -a.m_lb;
    c = a;
    for( int i = -10; i < 10; i++ )
    {
        arg[0] = i/10.0;
        for( int j = -10; j < 10; j++ )
        {
            arg[1] = j/10.0;
            // polynomial expansion
            double r = func.eval( arg );
            a.m_lb = min( a.m_lb, r );
            a.m_ub = max( a.m_ub, r );
            // actual function
            r = sin(arg[0]/2)/(2+cos(arg[1]/2+arg[0]*arg[0]));
            c.m_lb = min( c.m_lb, r );
            c.m_ub = max( c.m_ub, r );
        }
    }

    // DA bounding
    b = func.bound( );

    cout.precision( 16 );
    cout << "func:" << endl << func;
    cout << "Bounds:" << endl
         << "DA bound:       [" << b.m_lb << "," << b.m_ub << "]" << endl
         << "DA raster:      [" << a.m_lb << "," << a.m_ub << "]" << endl
         << "double raster:  [" << c.m_lb << "," << c.m_ub << "]" << endl << endl;
}


// Exercise 7.2.2: PID simulator with uncertain mass (DA)

// Model of controlled pendulum with uncertain mass
// x = ( theta, theta_dot, u )
AlgebraicVector<DA> pendulumRHSmass( AlgebraicVector<DA> x, double t )
{
    // pendulum constants
    static const double l = 0.61;       // m
    static const DA m = 0.21*(1+0.1*DA(1));       // kg
    static const double M = 0.4926;     // kg
    static const double g = 9.81;       // kg*m/s^2

    AlgebraicVector<DA> res(3);
    DA sint, cost;
    sint = sin(x[0]);
    cost = cos(x[0]);

    res[0] = x[1];
    res[1] = (x[2]+(M+m)*g*sint-m*l*sqr(x[1])*sint*cost)/((M+m)*l+m*l*sqr(cost));
    res[2] = 0;         // u is assumed constant unless changed externally by the controller

    return res;
}

void ex7_2_2( )
{
    const double Kp = 8, Ti = 3, Td = 0.3;                  // PID parameters
    const double setPt = 0.0, dt = 0.05;                    // Set point, controller time step (50 ms)
    const double umax = 10;                                 // maximum control (Exercise 7.2.1)

    ofstream file( "ex7_2_2.dat" );

    DA intErr = 0, lastErr = 0;                             // PID controller variables

    double t = 0;
    AlgebraicVector<DA> x(3);
    x[0] = 0.1; x[1] = 0; x[2] = 0;               // Initial condition

    // propagate the model state for 100 sec
    while( (t<40.0) && (abs(cons(x[0]))<1.5) )
    {
        // compute the PID control
        DA err = setPt-x[0];
        DA derr = (err-lastErr)/dt;
        intErr += lastErr*dt;
        lastErr = err;
        x[2] = Kp*(err + Td*derr + intErr/Ti);
        // prevent control saturation (Exercise 7.2.1)
//        x[2] = tanh(x[2]/umax)*umax;

        // output and propagate one time step (Exercise 7.1.4)
        Interval bx = x[0].bound( );
        Interval bu = x[2].bound( );
        file << t
        << "   " << cons(x[0]) << "   " << bx.m_ub << "   " << bx.m_lb
        << "   " << x[0].evalScalar(-1.0) << "   " << x[0].evalScalar(1.0)
        << "   " << cons(x[2]) << "   " << bu.m_ub << "   " << bu.m_lb
        << "   " << x[2].evalScalar(-1.0) << "   " << x[2].evalScalar(1.0) << endl;
        x = rk4( x, t, t+dt, pendulumRHSmass );
        t += dt;
    }

    cout << "Final angle theta:" << x << endl;
    if(abs(cons(x[0]))>1.5)
        cout << "WHOOPSY: Fell over after " << t << " seconds." << endl;

    file.close();
}



int main( int argc, char** argv )
{
    DA::init( 10, 2 );

    ex7_1_2( );
    ex7_1_3( );
    ex7_1_5( );

    ex7_2_2( );

    return 0;
}
