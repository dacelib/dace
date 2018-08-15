#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;


int main( void )
{
    
    DA::init( 10, 2 );
    
    const double tol = 1.0e-12;
    
    const double pi = 4.0*atan(1.0);
    const double mu = 1.0;
    
    DA a = 1.0;
    DA e = 0.5;
    double t = pi/2.0;
    
    DA M = sqrt(mu/(a*a*a))*t; // real at this stage
    
    DA EccAn = M; //first guess
    
    double err = abs(EccAn - e*sin(EccAn) - M);
    
    // Newton's method for the reference solution
    while(err>tol)
    {
        EccAn = EccAn - (EccAn - e*sin(EccAn) - M)/(1 - e*cos(EccAn));
        
        err = abs(EccAn - e*sin(EccAn) - M);
        
    }
    
    cout << setprecision(12);
    cout << "Reference solution: E = " << cons(EccAn) << endl;
    
    getchar();
    
    a = 1.0 + DA(1);
    e = 0.5 + DA(2);
    
    M = sqrt(mu/(a*a*a))*t; // now M is a DA
    
    // Newton's method for the Taylor expansion of the solution
    int i = 1;
    while( i <= 10){
        
        EccAn = EccAn - (EccAn - e*sin(EccAn) - M)/(1 - e*cos(EccAn));
        i *= 2;
        
    }
    
    cout << "Taylor expansion of E" << cons(EccAn) << endl;
    cout << EccAn << endl;
    
    getchar();
    
    cout << "Let's verify it is the Taylor expansion of the solution:" << endl;
    cout << "Evaluate (E - e*sin(E) - M) in DA" << endl;
    
    cout << EccAn - e*sin(EccAn) - M << endl;
    
}











