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
    
    DA y1 = sqr(sin(x)); // compute sin(x)^2
    
    cout << "sin(x)^2" << endl << y1 << endl;
    
    DA y2 = sqr(cos(x));  // compute cos(x)^2
    
    cout << "cos(x)^2" << endl << y2 << endl;
    
    // compute and print sin(x)^2+cos(x)^2
    cout << "sin(x)^2+cos(x)^2" << endl << y1+y2 << endl;
    
}