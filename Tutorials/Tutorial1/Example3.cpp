#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;


int main( void )
{
    
    DA::init( 1, 1 ); // Initialize DACE for 1st-order computations in 1 variable
    
    DA x = 3 + DA(1); // Initialize x as DA around 3
    
    cout << "x" << endl << x << endl;
    
    DA f = 1/(x+1/x); // Evaluate f(x) = 1/(x+1/x)
    
    cout << "f(x) = 1/(x+1/x)" << endl << f << endl;
    
}