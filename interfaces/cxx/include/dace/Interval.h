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
 * Interval.h
 *
 *  Created on: Mar 14, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_INTERVAL_H_
#define DINAMICA_INTERVAL_H_

#include "dace/Def.h"

namespace DACE{

/*! Class representing an interval. */
class DACE_API Interval
{
public:
    double m_lb;            //!< Lower bound.
    double m_ub;            //!< Upper bound.
};

}
#endif /* DINAMICA_INTERVAL_H_ */
