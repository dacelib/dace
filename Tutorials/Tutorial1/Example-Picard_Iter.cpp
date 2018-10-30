/// This is an example implementation of Picard iteration method for solving ODE

#include <dace/dace.h>
#include <cmath>

using namespace std;
using namespace DACE;

// The first example ODE
DA fun_1( DA x, DA t) {
    // dx = x + t, x(0) = 0
    // the analytic solution is: x = e^t - t - 1
    return x + t;
}

// The second example ODE
AlgebraicVector<DA> fun_2( AlgebraicVector<DA> x, DA t ) {
    // dx = y, dy = t
    // the analytic solution is:
    // x = x0 + y0*t + 1/6*t^3
    // y = y0 + 1/2*t^2
    AlgebraicVector<DA> dx;
    dx.resize(2);
    dx.at(0) = x.at(1);
    dx.at(1) = t; 

    return dx;
}

int main( void ) {
    DA::init( 12, 2 );

    // Example 1
    DA x0 = 0;
    DA t = DA(1);

    cout << "function: x = e^t - t - 1" << endl;
    DA xi, xj;
    xi = x0;
    for( int i = 0; i < 10; ++i ) {
        // iteration 10 times
        xj = x0 + fun_1(xi, t).integ(1);
        xi = xj;

        cout << "x_" << i+1 << endl << xj << endl;
    }

    // Example 2
    AlgebraicVector<DA> x_0(2);
    x_0.at(0) = 1;
    x_0.at(1) = 1;

    cout << "function: x = x0 + y0*t + 1/6*t^3, y = y0 + 1/2*t^2" << endl;
    AlgebraicVector<DA> x_i(2), x_j(2);
    x_i = x_0;
    for( int i = 0; i < 5; ++i ) {
        // iteration 5 times
        x_j = x_0 + fun_2(x_i, t).integ(1);
        x_i = x_j;

        cout << "iteration " << i+1 << endl;
        cout << "x" << endl << x_i.at(0) << endl;
        cout << "y" << endl << x_i.at(1) << endl;
    }

    return 0;
}