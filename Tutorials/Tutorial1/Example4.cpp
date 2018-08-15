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
    
    DA y = cos(x)-1; // Compute [cos(x)-1]
    
    cout << "[cos(x)-1]" << endl << y << endl;
    
    // Compute [cos(x)-1]^11
    for ( int i = 0; i < 10; i++)
    {
    	y = y*(cos(x)-1);
    }
    
    cout << "[cos(x)-1]^11" << endl << y << endl;
    
}