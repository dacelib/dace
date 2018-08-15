#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;


int main( void )
{
    
    DA::init( 10, 1 );
    
    DA x = DA(1);
    
    AlgebraicVector<DA> y(1), inv_y(1);
    
    y[0] = sin(x); // Compute Taylor expansion of sin(x)
    
    inv_y = y.invert(); // Invert Taylor polynomial
    
    // Compare with asin(x)
    cout << "Polynomial inversion of sin(x)" << endl << endl;
    cout << inv_y << endl << endl;
    
    cout << "asin(x)" << endl << endl;
    cout << asin(x) << endl;
    
}











