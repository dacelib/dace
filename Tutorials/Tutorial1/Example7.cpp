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
    
    DA::init( 10, 2 ); // initialize DACE for 10th-order computations in 2 variables
    
    cout << "Initialize x as two-dim DA vector around (0,0)" << endl << endl;
    
    AlgebraicVector<DA> x(2);
    
    x[0] = DA(1);
    x[1] = DA(2);
    
    cout << "x" << endl << x << endl;
    
    getchar();
    
    cout << "Evaluate sombrero function" << endl << endl;
    
    DA z = somb(x); // Evaluate sombrero function
    
    cout << "z" << endl << z << endl;
    
}











