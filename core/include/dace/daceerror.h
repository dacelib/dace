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
 *  daceerror.h
 *
 *  Created on: November 18, 2016
 *      Author: Politecnico di Milano
 */

#ifndef DINAMICA_DACEERROR_H
#define DINAMICA_DACEERROR_H

/** \addtogroup DACE Core
 *  @{
 */
typedef struct {
    int ID;     //!< Internal ID of the error
    const char* msg;    //!< Human readable error message
} errstrings;


DACE_API const errstrings DACEerr[] = {
    {   0, "Unknown DACE error. Contact DACE developers for filing a bug report."},
    {1001, "Dynamic memory allocation failure"},
    {1002, "Out of memory"},
    {1003, "DACE has not been initialized"},
    {1004, "DA object not allocated"},
    {1005, "Incorrect number of monomials"},
    {1006, "Incorrect DA coding arrays"},
    {1007, "Requested length too long"},
    {1008, "Error in monomial evaluation tree construction"},
    {   9, ""},
    {  10, ""},
    { 911, "Order and/or variable too large"},
    {  12, ""},
    {  13, ""},
    {  14, ""},
    {  15, ""},
    {  16, ""},
    {  17, ""},
    {  18, ""},
    {  19, ""},
    {  20, ""},
    { 621, "Not enough storage"},
    { 622, "Order too large"},
    { 623, "Number of variables too high"},
    { 624, "Invalid independent variable"},
    { 625, "Invalid DA codes"},
    { 626, "Invalid encoded exponent"},
    {  27, ""},
    {  28, ""},
    {  29, ""},
    {  30, ""},
    { 631, "Invalid data"},
    { 632, "Unknown format"},
    { 633, "DA vector too long"},
    { 634, "Not enough lines to read"},
    {  35, ""},
    {  36, ""},
    {  37, ""},
    {  38, ""},
    {  39, ""},
    {  40, ""},
    { 641, "Dividing by zero"},
    { 642, "Inverse does not exists"},
    { 643, "Non-integer power of non-positive DA"},
    { 644, "Zero-th root does not exist"},
    { 645, "Even root of negative DA"},
    { 646, "Odd root of zero DA"},
    { 647, "Negative constant part in logarithm"},
    { 648, "Base of logarithm must be positive"},
    { 649, "Cosine is zero in tangent"},
    { 650, "Out of domain"},
    { 651, "No estimate is possible"},
    {  52, ""},
    {  53, ""},
    {  54, ""},
    {  55, ""},
    {  56, ""},
    {  57, ""},
    {  58, ""},
    {  59, ""},
    {  60, ""},
    { 161, "Free or invalid variable"},
    { 162, "Truncation order too high"},
    { 163, "Inacurate estimate"},
    { 164, "Numbering out of order"},
    { 165, "Too many variables"},
    { 166, "Duplicate monomial"},
    { 167, "Order increased to 1"},
    { 168, "Variable increased to 1"},
    {  69, ""},
    {  70, ""}
}; //!< Variable containing all errors strings and codes
/** @}*/
#endif
