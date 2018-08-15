#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;


template<typename T> T ErrFunc(T x)
{
	const double pi = 4.0*atan(1.0);
	
	T z = 1.0/sqrt(2.0*pi)*exp(-x*x/2);
	return z;
}


int main( void )
{
    
    DA::init( 24, 1 ); // initialize DACE for 20th-order computations in 1 variable
    
    DA x = DA(1);
    
    DA y = ErrFunc(x); // compute Taylor expansion of the erf
    
    DA Inty = y.integ(1); // compute the Taylor expansion of the indefinite integral of erf
    
    double value = Inty.evalScalar(1.0)-Inty.evalScalar(-1.0); // compute int_{-1}^{+1} (erf)
    
    cout << "int_{-1}^{+1} (erf)" << endl;
    cout << setprecision(12);
    cout << "Exact result: " << 0.682689492137 << endl;
    cout << "Approx result " << value << endl;
    
}