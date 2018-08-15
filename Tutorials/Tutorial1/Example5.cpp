#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;


int main( void )
{
    
    DA::init( 20, 1 );
    
    DA x = DA(1);
    
    DA y = sin(x);
    
    DA dy = y.deriv(1); // compute Taylor expansion of d[sin(x)]/dx
    
    // print d[sin(x)]/dx and cos(x) to compare
    cout << "d[sin(x)]/dx" << endl << dy << endl;
    cout << "cos(x)" << endl << cos(x) << endl;
    
    getchar();
    
    DA int_y = y.integ(1); // compute Taylor expansion of int[sin(x)dx]
    
    // print int[sin(x)dx] and -cos(x) to compare
    cout << "int[sin(x)dx]" << endl << int_y << endl;
    cout << "-cos(x)" << endl << -cos(x) << endl;
    
}