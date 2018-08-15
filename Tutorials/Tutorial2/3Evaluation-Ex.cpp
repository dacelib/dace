#include <dace/dace.h>
#include <cmath>
#include <fstream>
using namespace std; using namespace DACE;

const double pi =  3.141592653589793;
const double val1 = 1.685401585899429; // reference value of the integral over -1, +1
const double val2 = 1.990644530037905; // reference value of the integral over -2, +2

// Exercise 3.1.1: plot a 1D polynomial
void ex3_1_1( )
{
    const double x0 = 0.0;          // expansion point
    const int N = 100;              // number of points in grid
    const double hw = 2.0;          // length of grid in each direction from expansion point

    DA x = DA(1)+x0;
    DA func = exp(-x*x);
    ofstream file;

    file.open( "ex3_1_1.dat" );
    for( int i = 0; i < N; i++ )
    {
        double xx = -hw + i*2.0*hw/(N-1);   // point on the grid on [-hw,hw]
        double rda = func.evalScalar( xx ); // note: this is not an efficient way to repeatedly evaluate the same polynomial
        xx += x0;   // add expansion point x0 for double evaluation
        double rdouble = exp(-xx*xx);

        file << xx << "   " << rda << "   " << rdouble << endl;
    }
    file.close( );
    // gnuplot command: plot 'ex3_1_1.dat'u 1:2 w l, 'ex3_1_1.dat' u 1:3 w l
    // or for the error: plot 'ex3_1_1.dat' u 1:($2-$3) w l
}

// Exercise 3.1.2: plot a 2D polynomial
template<typename T> T somb( T x, T y )
{
    T r = sqrt( x*x + y*y );
    return sin(r)/r;
}

void ex3_1_2( )
{
    const double x0 = 1.0;          // expansion point x
    const double y0 = 1.0;          // expansion point y
    const int N = 50;               // number of points in grid

    DA x = DA(1)+x0; DA y = DA(2)+y0;
    DA func = somb(x,y);
    vector<double> arg(2);          // vector holding two doubles
    ofstream file;

    file.open( "ex3_1_2.dat" );
    for( int i = 0; i < N; i++ )
    {
        arg[0] = -1.0 + i*2.0/(N-1);  // x coordinate on the grid on [-1,1]
        for( int j = 0; j < N; j++ )
        {
            arg[1] = -1.0 + j*2.0/(N-1);    // y coordinate on the grid on [-1,1]
            double rda = func.eval( arg );  // note: this is not an efficient way to repeatedly evaluate the same polynomial
            double rdouble = somb( x0+arg[0], y0+arg[1] );
            file << arg[0] << "   " << arg[1] << "   " << rda << "   " << rdouble << endl;
        }
        file << endl;   // empty line between lines of data for gnuplot
    }
    file.close( );
    // gnuplot command: splot 'ex3_1_2.dat' u 1:2:3 w l, 'ex3_1_2.dat' u 1:2:4 w l
    // or for the error: splot 'ex3_1_2.dat' u 1:2:($3-$4) w l
}

// Exercise 3.1.3: Sinusitis
void ex3_1_3( )
{
    DA::pushTO( 10 );

    DA x = DA(1);
    DA sinda = sin(x);
    DA res1 = sin( x+2 );                 // compute directly sin(1+DA(1))
    DA res2 = sinda.evalScalar( x+2 );    // evaluate expansion of sine with 1+DA(1)

    cout << "Exercise 4.1.3: Sinusitis" << endl << res1-res2 << endl;

    DA::popTO( );
}

// Exercise 3.1.4: Gauss integral I
void ex3_1_4( void )
{
    ofstream file;

    file.open( "ex3_1_4.dat" );
    file.precision( 16 );           // print 16 digits, similar to MATLAB's "format long"
    for( int order = 1; order <= 40; order++ )
    {
        DA::setTO( order );                             // limit the computation order
        DA t = DA(1);
        DA erf = 2.0/sqrt(pi)*exp(-sqr(t)).integ(1);    // error function erf(x)
        double res = erf.evalScalar( 1.0 )-erf.evalScalar( -1.0 );
        file << order << "   " << res << "   " << log(abs(res-val1))/log(10.0) << endl;
    }
    file.close( );
    // gnuplot command: plot 'ex3_1_4.dat'u 1:2 w l
    // or for the error: plot 'ex3_1_4.dat'u 1:3 w l
}

// Exercise 3.2.1 & 3.2.2: Gauss integral II
double gaussInt( double a, double b )       // compute integral of Gaussian on interval [a,b]
{
    DA t = (a+b)/2.0+DA(1);     // expand around center point
    DA func = 2.0/sqrt(pi)*exp(-sqr(t)).integ(1);
    return func.evalScalar( (b-a)/2.0 ) - func.evalScalar( -(b-a)/2.0 );    // evaluate over -+ half width
}

void ex3_2_1( void )
{
    const double hw = 2.0;      // half-width of the interval to integrate on, i.e. [-hw,hw]
    ofstream file;
    double res;

    file.open( "ex3_2_1.dat" );
    file.precision( 16 );
    DA::pushTO( 9 );
    for( int N = 1; N <= 30; N++ )
    {
        res = 0.0;
        for( int i = 1; i <= N; i++ )
        {
            double ai  = -hw + (i-1)*2.0*hw/N;
            double ai1 = -hw + i*2.0*hw/N;
            res += gaussInt( ai, ai1 );
        }
        file << N << "   " << res << "   " << log(abs(res-val2))/log(10.0) << endl;
    }
    DA::popTO( );
    // compare to single expansion at full computation order
    res = gaussInt( -hw, hw );
    file << endl << 1 << "   " << res << "   " << log(abs(res-val2))/log(10.0) << endl;
    file.close( );
    // gnuplot command: plot 'ex3_2_1.dat'u 1:3 w lp
}

int main( void )
{
    DA::init( 30, 2 );      // init with maximum computation order

    ex3_1_1( );
    ex3_1_2( );
    ex3_1_3( );
    ex3_1_4( );

    ex3_2_1( );
}

