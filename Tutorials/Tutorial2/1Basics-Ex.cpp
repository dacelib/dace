#include <dace/dace.h>
#include <cmath>

using namespace std;
using namespace DACE;

// Exercise 1.1.2: first steps
void ex1_1_2( )
{
    DA x = DA(1);
    DA func = 3*(x+3)-x+1-(x+8);
    cout << "Exercise 1.1.2: First steps\n" << func << endl;
}

// Exercise 1.1.3: different expansion point
void ex1_1_3( )
{
    DA x = DA(1);
    DA func = sin(1.0+x);
    cout << "Exercise 1.1.3: Different expansion point\n" << func << endl;
}

// Exercise 1.1.4: a higer power
void ex1_1_4( )
{
    DA x = DA(1);
    DA func = sin(x);
    DA res = 1.0;       // this makes res a constant function P(x) = 1.0
    for( int i = 0; i < 11; i++ ) res = res*sin(x);
    cout << "Exercise 1.1.4: A higher power\n" << res << endl;
}

// Exercise 1.1.5: two arguments
template<typename T> T ex1_1_5( T x, T y )
{
    return sqrt( 1.0 + x*x + y*y );
}

// Exercise 1.2.1: identity crisis
void ex1_2_1( )
{
    DA x = DA(1);
    DA s2 = sin(x)*sin(x);
    DA c2 = cos(x)*cos(x);
    cout << "Exercise 1.2.1: Identity crisis\n" << s2 << c2 << s2+c2 << endl;
}

// Exercise 1.2.2: Breaking bad
template<typename T> T ex1_2_2( T x, T y )
{
    T r = sqrt( x*x + y*y );
    return sin(r)/r;
}

int main( void )
{
    DA::init( 10, 2 );
    cout.precision( 15 );

    DA x = DA(1);
    DA y = DA(2);

    ex1_1_2( );
    ex1_1_3( );
    ex1_1_4( );
    cout << "Exercise 1.1.5: Two arguments\n" << ex1_1_5( x, y ) << endl;

    ex1_2_1( );
    try {
    cout << "Exercise 1.2.2: Breaking bad\n"
         << ex1_2_2( 0.0, 0.0 ) << endl
         << ex1_2_2( x, y ) << endl;
    } catch( exception const& ex ) { cout << ex.what( ) << endl; }
}