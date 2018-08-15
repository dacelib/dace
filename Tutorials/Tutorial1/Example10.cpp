#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;


int main( void )
{
    
    DA::init( 10, 3 );
    
    const double pi = 4.0*atan(1.0);
    
    // initialize cyl
    AlgebraicVector<DA> cyl(3);
    
    cyl[0] = 100.0+DA(1);
    cyl[1] = (0.0+DA(2))*pi/180.0;
    cyl[2] = 0.0 + DA(3);
    
    // initialize cart and compute transformation
    AlgebraicVector<DA> cart(3);
    
    cart[0] = cyl[0]*cos(cyl[1]);
    cart[1] = cyl[0]*sin(cyl[1]);;
    cart[2] = cyl[2];
    
    // subtract constant part to build DirMap
    AlgebraicVector<DA> DirMap(3);
    
    DirMap = cart-cart.cons();
    
    cout << "Direct map: from cylindric to cartesian (DirMap)" << endl << endl;
    cout << DirMap << endl << endl;
    
    // Invert DirMap to obtain InvMap
    AlgebraicVector<DA> InvMap(3);
    
    InvMap = DirMap.invert();
    
    cout << "Invers map: from cartesian to cylindric (InvMap)" << endl << endl;
    cout << InvMap << endl << endl;

    getchar();
    
    // Verification
    cout << "Concatenate the direct and inverse map: (DirMap) o (InvMap) = DirMap(InvMap) = I" << endl << endl;
    cout << DirMap.eval(InvMap) << endl;
    
}











