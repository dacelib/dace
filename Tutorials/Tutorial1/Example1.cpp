#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;

int main( void )
{
    
    DA::init( 20, 1 ); // initialize DACE for 20th-order computations in 1 variable
    
    DA x = DA(1); // initialize x as DA
    
    DA y = sin(x); // compute y = sin(x)
    
    // print x and y to screen
    cout << "x" << endl << x << endl;
    cout << "y = sin(x)" << endl << y;
    
}