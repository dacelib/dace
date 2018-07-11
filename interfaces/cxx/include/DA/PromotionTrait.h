/******************************************************************************
*                                                                             *
* DIFFERENTIAL ALGEBRA CORE ENGINE                                            *
*                                                                             *
*******************************************************************************
*                                                                             *
* Copyright 2016 Politecnico di Milano (2014 Dinamica Srl)                    *
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
 * PromotionTrait.h
 *
 *  Created on: Sep. 15, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_PROMOTIONTRAIT_H_
#define DINAMICA_PROMOTIONTRAIT_H_

namespace DACE{

// forward declaration
class DA;

/***********************************************************************************
      PromotionTrait template class
************************************************************************************/
/* General template implementation:
   By default, everything is a double and two things of the same type return that type. */
template<class T1, class T2> class PromotionTrait { public: typedef double returnType; };
template<class T> class PromotionTrait<T,T> { public: typedef T returnType; };

// type A with B (and B with A) yields a C
#define ADD_PROMOTION(A,B,C) \
    template<> class PromotionTrait<A,B> { public: typedef C returnType; }; \
    template<> class PromotionTrait<B,A> { public: typedef C returnType; };

// type A with any other type (and any other type with A) yields a C
// includes also the case A with A to resolve the problem of which template to use otherwise
#define ADD_PROMOTION_ALL(A,C) \
    template<typename T> class PromotionTrait<A,T> { public: typedef C returnType; }; \
    template<typename T> class PromotionTrait<T,A> { public: typedef C returnType; }; \
    template<> class PromotionTrait<A,A> { public: typedef C returnType; };

/* Specialization for DA: everything with DA returns a DA */
ADD_PROMOTION_ALL(DA,DA)

#undef ADD_PROMOTION
#undef ADD_PROMOTION_ALL
}

// not sure this is a good idea
//#define PROMOTE(A,B) typename PromotionTrait<A,B>::returnType

#endif /* DINAMICA_PROMOTIONTRAIT_H_ */
