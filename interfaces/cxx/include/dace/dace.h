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
 * dace.h
 *
 *  Created on: Jan 14, 2015
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DACE_H_
#define DINAMICA_DACE_H_

// This file just brings in all the headers of the DACE related classes and functions
#include "dace/config.h"
#include "dace/PromotionTrait.h"
#include "dace/MathExtension.h"
#include "dace/DACEException.h"
#include "dace/Monomial.h"
#include "dace/Interval.h"
#include "dace/DAFormatter.h"
#include "dace/compiledDA.h"
#include "dace/DA.h"
#include "dace/AlgebraicVector.h"
#ifdef WITH_ALGEBRAICMATRIX
#include "dace/AlgebraicMatrix.h"
#endif /* WITH_ALGEBRAICMATRIX */

// include the template implementations here at the end after everything is properly defined
#include "dace/compiledDA_t.h"
#include "dace/DA_t.h"
#include "dace/AlgebraicVector_t.h"
#ifdef WITH_ALGEBRAICMATRIX
#include "dace/AlgebraicMatrix_t.h"
#endif /* WITH_ALGEBRAICMATRIX */

#endif /* DINAMICA_DACE_H_ */
