/******************************************************************************
*                                                                             *
* DIFFERENTIAL ALGEBRA CORE ENGINE                                            *
*                                                                             *
*******************************************************************************
*                                                                             *
* Copyright 2025 Alexander Wittig                                             *
* Licensed under the Apache License, Version 2.0 (the "License");             *
* you may not use this file except in compliance with the License.            *
* You may obtain a copy of the License at                                     *
*                                                                             *
*    http://www.apache.org/licenses/LICENSE-2.0                               *
*                                                                             *
* Unless required by applicable law or agreed to in writing, software         *
* distributed under the License is distributed on an "AS IS" BASIS,           *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    *
* See the License for the specific language governing permissions and         *
* limitations under the License.                                              *
*                                                                             *
*******************************************************************************/

/*
 * compat_boost_odeint.h
 *
 *  Created on: Apr 25, 2025
 *      Author: Alexander Wittig
 */

#ifndef DINAMICA_COMPAT_BOOST_ODEINT_H_
#define DINAMICA_COMPAT_BOOST_ODEINT_H_

/*! This file needs to be included after boost/numeric/odeint.hpp and dace/dace.h to provide
    a compatibility shim to allow odeint to work with DACE::AlgebraicVector<> as a state type.
    It can be used both with DA and double data types.

    Additionally, after including this file, an implementation of abs(DACE::DA) must be provided
    within the DACE namespace.
    For most practical uses, this can be done by selecting one of the predefined implementations
    provided in DACE::abs_cons, DACE::abs_max, or DACE::abs_sum.

    Example:
    \code
        #include <boost/numeric/odeint.hpp>
        #include <dace/dace.h>
        #include <dace/compat_boost_odeint.h>

        // select implementation of abs(DA)
        namespace DACE { using DACE::abs_max::abs; }


        using namespace boost::numeric::odeint;
        using namespace DACE;

        // define RHS, call DA::init, ...

        typedef state_type AlgebraicVector<DA>;

        state_type x(3);
        x[0] = 1.0; x[1] = 2.0; x[3] = 3.0;
        x += 0.01 * x.identity(3);
        integrate_adaptive( make_controlled( 1e-8, 1e-8, runge_kutta_dopri5<state_type>() ), RHS, x, 0.0, 10.0, 0.1 );
    \endcode
 */

 namespace boost { namespace numeric { namespace odeint {
    // mark AlgebraicVectors as resizable (using the standard container interface inherited from std::vector)
    template<typename T> struct is_resizeable<DACE::AlgebraicVector<T>>
    {
        typedef boost::true_type type;
        static const bool value = type::value;
    };
} } }

#endif /* DINAMICA_COMPAT_BOOST_ODEINT_H_ */
