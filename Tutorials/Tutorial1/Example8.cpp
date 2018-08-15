#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;


template<typename T> T somb(AlgebraicVector<T> x)
{
    
	return sin(sqrt(sqr(x[0]) + sqr(x[1])))/sqrt(sqr(x[0]) + sqr(x[1]));

}


int main( void )
{
    
    DA::init( 1, 2 ); // initialize DACE for 1st-order computations in 2 variables
    
    cout << "Initialize x as two-dim DA vector around (2,3)" << endl << endl;
    
    AlgebraicVector<DA> x(2);
    
    x[0] = 2.0 + DA(1);
    x[1] = 3.0 + DA(2);
    
    DA z = somb(x); // Evaluate sombrero function
    
    cout << "Sombrero function" << endl << endl;
    cout << z << endl << endl;
    
    AlgebraicVector<DA> grad_z(2); // Declare DA vector that will contain the gradient of sombrero function
    
    grad_z = z.gradient(); // Compute gradient of sombrero function
    
    cout << "Gradient of sombrero function" << endl << endl;
    cout << grad_z << endl;
    
}











