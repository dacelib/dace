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
 * AlgebraicMatrix_t.h
 *
 *  Created on: July 17, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DAMATRIX_T_H_
#define DINAMICA_DAMATRIX_T_H_

// C++ stdlib classes used only internally in this implementation
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <string>

// DACE classes
#include "dace/PromotionTrait.h"
#include "dace/AlgebraicVector.h"
#include "dace/AlgebraicMatrix.h"

namespace DACE{

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicVector<U> &obj1, const AlgebraicMatrix<V> &obj2 ) {
/*! Compute the multiplication between a AlgebraicVector (row) and a AlgebraicMatrix.
   \param[in] obj1 a AlgebraicVector
   \param[in] obj2 a AlgebraicMatrix
   \return A new AlgebraicVector, containing the result of the operation (obj1*obj2)
   \sa AlgebraicMatrix<T>::operator*
 */
    // Check that vector length and number of rows of matrix are equal
    if (obj1.size()!=obj2.nrows())
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::operator*: objects in vector/matrix multiplication must have compatible size.");

    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(obj2.ncols());

    for(unsigned int i=0; i<obj2.ncols(); i++){
        temp[i] = 0.;
        for (unsigned int j=0; j<obj2.nrows();j++)
        temp[i] += obj1[j]*obj2.at(j,i);
    }

    return temp;
}

template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicMatrix<U> &obj1, const AlgebraicVector<V> &obj2 ){
/*! Compute the multiplication between a AlgebraicMatrix and a AlgebraicVector (column) .
   \param[in] obj1 a AlgebraicMatrix.
   \param[in] obj2 a AlgebraicVector.
   \return A new AlgebraicVector, containing the result of the operation (obj1*obj2).
   \sa AlgebraicMatrix<T>::operator*
 */
    // Check that vector length and number of columns of matrix are equal
    if (obj1.ncols()!=obj2.size())
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::operator*: objects in vector/matrix multiplication must have compatible size.");

    AlgebraicVector<typename PromotionTrait< U, V >::returnType> temp(obj1.nrows());

    for(unsigned int i=0; i<obj1.nrows(); i++){
        temp[i] = 0.;
        for (unsigned int j=0; j<obj1.ncols();j++)
            temp[i] += obj1.at(i,j)*obj2[j];
    }

    return temp;
}

template<class T> void AlgebraicMatrix<T>::resize(int size) {
/*! Resize AlgeraicMatrix to a square AlgebraicMatrix of size.
   The original values are kept at the original location if they are inside bounds of the new matrix.
   \param[in] size Size of the matrix (number of rows/columns).
   \sa AlgebraicMatrix<T>::resize
 */
    // store initial data
    unsigned int oldrows = this->nrows();
    unsigned int oldcols = this->ncols();

    std::vector<T> temp(this->_data);

    this->_nrows = size;
    this->_ncols = size;
    this->_data.resize(size*size);

    // Sort old values to keep same position of original matrix
    for(unsigned int irow=0; irow<(unsigned)size; irow++)
        for(unsigned int icol=0; icol<(unsigned)size; icol++)
            if ((irow<oldrows)&&(icol<oldcols))
                this->_data[irow*this->_ncols + icol] = temp[irow*oldcols+icol];
            else
                this->_data[irow*this->_ncols + icol] = 0.;
}

template<class T> void AlgebraicMatrix<T>::resize(int rows, int cols) {
/*! Resize AlgeraicMatrix to a rectangular AlgebraicMatrix with size rows\f$\times\f$ cols.
   The original values are kept at the original location if they are inside bounds of the new matrix.
   \param[in] rows Number of rows of the resized AlgebraicMatrix.
   \param[in] cols Number of columns of the resized AlgebraicMatrix.
   \sa AlgebraicMatrix<T>::resize
 */
    // store initial data
    unsigned int oldrows = this->nrows();
    unsigned int oldcols = this->ncols();
    std::vector<T> temp(this->_data);

     this->_nrows = rows;
     this->_ncols = cols;
     this->_data.resize(rows * cols);

    // Sort old values to keep same position of original matrix
    for(unsigned int irow=0; irow<(unsigned)rows; irow++)
        for(unsigned int icol=0; icol<(unsigned)cols; icol++)
            if ((irow<oldrows)&&(icol<oldcols))
                this->_data[irow*cols + icol] = temp[irow*oldcols+icol];
            else
                this->_data[irow*cols + icol] = 0.;
}

/***********************************************************************************
*     Element access routines
************************************************************************************/
template<class T> T& AlgebraicMatrix<T>::at(const unsigned int irow, const unsigned int icol) {
/*! Reads/write element from/to AlgebraicMatrix.
   \param[in] irow row index.
   \param[in] icol column index.
   \return The element of the AlgebraicMatrix.
   \sa AlgebraicMatrix<T>::at
 */
    if ( !(irow<(this->_nrows) )&&!( icol<(this->_ncols) ))
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::at: matrix element position out of bound.");

    return this->_data[irow*this->_ncols + icol];
}

template<class T> const T& AlgebraicMatrix<T>::at(const unsigned int irow, const unsigned int icol) const {
/*! Reads/write element from/to AlgebraicMatrix.
   \param[in] irow row index.
   \param[in] icol column index.
   \return The element of the AlgebraicMatrix.
   \sa AlgebraicMatrix<T>::at
 */
    if ( !(irow<(this->_nrows) )&&!( icol<(this->_ncols)) )
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::at: matrix element position out of bound.");

    return this->_data[irow*this->_ncols + icol];
}

template<class T> std::vector<T> AlgebraicMatrix<T>::getrow(const unsigned int irow) const {
/*! Extracts a row of an AlgebraicMatrix.
    The result is copied in a new std::vector.
   \param[in] irow index of the row to be extracted
   \return A new std::vector containing the row of the AlgebraicMatrix.
   \sa AlgebraicMatrix<T>::getrow
 */
    // check that desired row is within matrix bounds
    if( !(irow<(this->_nrows)) )
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::row: row out of bound.");

    // Declare a vector with length equal to number of columns of the matrix
    std::vector<T> temp(this->_ncols);

    // Extract row and store elements into vector
    for(unsigned int i=0; i<(this->_ncols); i++)
    {
        temp[i] = this->_data[irow*this->_ncols+i];
    }

    return temp;
}

template<class T> std::vector<T> AlgebraicMatrix<T>::getcol(const unsigned int icol) const {
/*! Extracts a column of an AlgebraicMatrix.
    The result is copied in a new std::vector.
   \param[in] icol index of the column to be extracted
   \return A new std::vector containing the column of the AlgebraicMatrix.
   \sa AlgebraicMatrix<T>::getcol
 */
    // check that desired column is within matrix bounds
    if( !(icol<(this->_ncols)) )
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::col: column out of bound.");

    // Declare a vector with length equal to number of rows of the matrix
    std::vector<T> temp(this->_nrows);

    // Extract column and store elements into vector
    for(unsigned int i=0; i<(this->_nrows); i++)
    {
        temp[i] = this->_data[i*this->_ncols+icol];
    }

    return temp;
}

template<class T> void AlgebraicMatrix<T>::setrow(const unsigned int irow, const std::vector< T >& obj) {
/*! Insert std::vector into row of AlgebraicMatrix of the same type.
   \param[in] irow row to be written
   \param[in] obj  std::vector to be inserted as a row
   \sa AlgebraicMatrix<T>::setrow
 */
    // check that length of vector and number of columns match
    if (obj.size()!=this->_ncols)
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::setrow: vector too large to be stored in matrix.");

    // Insert elements at desired row
    for(unsigned int j=0; j<(this->_ncols); j++)
        this->_data[irow*this->_ncols+j] = obj[j];
}

template<class T>  void AlgebraicMatrix<T>::setcol(const unsigned int icol, const std::vector< T >& obj) {
/*! Insert std::vector into column of AlgebraicMatrix of the same type.
   \param[in] icol column to be written
   \param[in] obj  std::vector to be written
   \sa AlgebraicMatrix<T>::setrow
 */
    // check that length of vector and number of rows match
    if (obj.size()!=this->_nrows)
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::setcol: vector too large to be stored in matrix.");

    // Insert elements at desired column
    for(unsigned int i=0; i<(this->_nrows); i++)
        this->_data[i*this->_ncols+icol] = obj[i];
}

template<class T> AlgebraicMatrix<T> AlgebraicMatrix<T>::submat(const unsigned int first_row, const unsigned int first_col, const unsigned int last_row, const unsigned int last_col) const {
/*! Extracts submatrix of AlgebraicMatrix. The result is stored into a new AlgebraicMatrix.
   \param[in] first_row index of the first row to be extracted.
   \param[in] last_row  index of the last row to be extracted.
   \param[in] first_col index of the first column to be extracted.
   \param[in] last_col  index of the last column to be extracted.
   \return AlgebraicMatrix containing the desired submatrix.
   \sa AlgebraicMatrix<T>::submat
 */
  // Check that desired last_row and last_col are actually larger than first_row and first_col
  if( (last_row<first_row)||(last_col<first_col) )
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::submat: first row/column index larger than last row/column index.");

  // Check that desired submatrix can be extracted from original matrix
  if( !(last_row<(this->_nrows))||!(last_col<(this->_ncols)) )
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::submat: last row/column index exceeds number of rows/columns.");

  unsigned int nrows = last_row-first_row+1;
  unsigned int ncols = last_col-first_col+1;

  // Declare temporary matrix for output
  AlgebraicMatrix<T> temp(nrows, ncols);

  for (unsigned int i=0; i<nrows; i++)
        for (unsigned int j=0; j<ncols; j++)
            temp.at(i,j) = this->at(i+first_row,j+first_col);

  return temp;
}

template<class T> AlgebraicMatrix<T> AlgebraicMatrix<T>::submat(const unsigned int last_row, const unsigned int last_col) const {
/*! Extracts submatrix of AlgebraicMatrix, starting from element (0,0).
    The result is stored into a new AlgebraicMatrix.
   \param[in] last_row  index of the last row to be extracted.
   \param[in] last_col  index of the last column to be extracted.
   \return AlgebraicMatrix containing the desired submatrix.
   \sa AlgebraicMatrix<T>::submat
 */
  return this->submat(0,0,last_row,last_col);
}

/***********************************************************************************
*     Algebraic operations
************************************************************************************/
template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator+( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2) {
/*! Compute the addition between two AlgebraicMatrices.
   \param[in] obj1 the first AlgebraicMatrix.
   \param[in] obj2 the second AlgebraicMatrix.
   \return A new AlgebraicMatrix, containing the result of the operation (obj1+obj2).
   \sa AlgebraicMatrix<T>::operator+
 */
    // check that the two matrices have the same size
    if((obj1.ncols()!=obj2.ncols())||(obj1.nrows()!=obj2.nrows()))
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::operator+: Inputs must have the same size, unless one is a scalar value.");

    unsigned int n_rows = obj1.nrows();
    unsigned int n_cols = obj1.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
            temp.at(i,j) = obj1.at(i,j) + obj2.at(i,j);

    return temp;
}

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator+( const AlgebraicMatrix<U> &obj1, const V &obj2 ) {
/*! Compute the addition between a AlgebraicMatrix and a scalar value.
   \param[in] obj1 a AlgebraicMatrix.
   \param[in] obj2 a scalar value.
   \return A new AlgebraicMatrix, containing the result of the operation (obj1+obj2).
   \sa AlgebraicMatrix<T>::operator+
 */
    unsigned int n_rows = obj1.nrows();
    unsigned int n_cols = obj1.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
      for (unsigned int j=0; j<n_cols; j++)
    temp.at(i,j) = obj1.at(i,j) + obj2;

    return temp;
}

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator+( const U &obj1, const AlgebraicMatrix<V> &obj2 ) {
/*! Compute the addition between a scalar value and a AlgebraicMatrix.
   \param[in] obj1 a scalar value.
   \param[in] obj2 a AlgebraicMatrix.
   \return A new AlgebraicMatrix, containing the result of the operation (obj1+obj2).
   \sa AlgebraicMatrix<T>::operator+
 */
    unsigned int n_rows = obj2.nrows();
    unsigned int n_cols = obj2.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
            temp.at(i,j) = obj1 + obj2.at(i,j);

    return temp;
}

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator-( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2) {
/*! Compute the subtraction between two AlgebraicMatrices.
   \param[in] obj1 first AlgebraicMatrx.
   \param[in] obj2 second AlgebraicMatrix.
   \return A new AlgebraicMatrix, containing the result of the operation (obj1-obj2).
   \sa AlgebraicMatrix<T>::operator-
 */
    // check that the two matrices have the same size
    if((obj1.ncols()!=obj2.ncols())||(obj1.nrows()!=obj2.nrows()))
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::operator-: Inputs must have the same size, unless one is a scalar value.");

    unsigned int n_rows = obj1.nrows();
    unsigned int n_cols = obj1.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
            temp.at(i,j) = obj1.at(i,j)-obj2.at(i,j);

    return temp;
}

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator-( const AlgebraicMatrix<U> &obj1, const V &obj2 ) {
/*! Compute the subtraction between a AlgebraicMatrix and a scalar value.
   \param[in] obj1 a AlgebraicMatrix.
   \param[in] obj2 a scalar value.
   \return A new AlgebraicMatrix, containing the result of the operation (obj1-obj2).
   \sa AlgebraicMatrix<T>::operator-
 */
    unsigned int n_rows = obj1.nrows();
    unsigned int n_cols = obj1.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
            temp.at(i,j) = obj1.at(i,j) - obj2;

    return temp;
}


template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator-( const U &obj1, const AlgebraicMatrix<V> &obj2 ) {
/*! Compute the subtraction between a scalar value and a AlgebraicMatrix.
   \param[in] obj1 a scalar value.
   \param[in] obj2 a AlgebraicMatrix.
   \return A new AlgebraicMatrix, containing the result of the operation (obj1-obj2).
   \sa AlgebraicMatrix<T>::operator-
 */
    unsigned int n_rows = obj2.nrows();
    unsigned int n_cols = obj2.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
            temp.at(i,j) = obj1 - obj2.at(i,j);

    return temp;
}

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2 ){
/*! Compute the multiplication between two AlgebraicMatrices.
* \param[in] obj1 first AlgebraicMatrix  \f$ (n \times m) \f$
* \param[in] obj2 second AlgebraicMatrix \f$ (m \times p) \f$.
* \return A new AlgebraicMatrix, of size \f$ (n \times p) \f$ containing the result of the operation (obj1*obj2).
* \sa AlgebraicMatrix<T>::operator*
*/
    // check that first matrix columns and second matrix row are equal
    if( obj1.ncols()!=obj2.nrows() )
    throw std::runtime_error("DACE::AlgebraicMatrix<T>::operator*: Number of columns of first matrix must be equal to num,ber of rows of second matrix.");

    unsigned int N = obj1.nrows();
    unsigned int M = obj1.ncols();
    unsigned int P = obj2.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(N,P);

    for (unsigned int i=0; i<N; i++)
        for (unsigned int j=0; j<P; j++) {
            temp.at(i,j) = 0;
            for (unsigned int k=0; k<M; k++)
                temp.at(i,j) += obj1.at(i,k)*obj2.at(k,j); }

    return temp;
}

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicMatrix<U> &obj1, const V &obj2 ) {
/*! Compute the multiplication between a AlgebraicMatrix and a scalar value.
   \param[in] obj1 a AlgebraicMatrix.
   \param[in] obj2 a scalar value.
   \return A new AlgebraicMatrix, containing the result of the operation (obj1*obj2).
   \sa AlgebraicMatrix<T>::operator*
 */
    unsigned int n_rows = obj1.nrows();
    unsigned int n_cols = obj1.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
        for (unsigned int j=0; j<n_cols; j++)
            temp.at(i,j) = obj1.at(i,j)*obj2;

    return temp;
}

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator*( const U &obj1, const AlgebraicMatrix<V> &obj2 ) {
/*! Compute the multiplication between a scalar value and a AlgebraicMatrix.
   \param[in] obj1 a scalar value.
   \param[in] obj2 a AlgebraicMatrix.
   \return A new DAvector, containing the result of the operation (obj1*obj2).
   \sa AlgebraicMatrix<T>::operator*
 */
    unsigned int n_rows = obj2.nrows();
    unsigned int n_cols = obj2.ncols();

    AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> temp(n_rows, n_cols);

    for (unsigned int i=0; i<n_rows; i++)
      for (unsigned int j=0; j<n_cols; j++)
    temp.at(i,j) = obj1*obj2.at(i,j);

    return temp;
}

/***********************************************************************************
*     Matrix operations
************************************************************************************/
template<class T> AlgebraicMatrix<T> AlgebraicMatrix<T>::transpose() const {
/*! Transpose matrix. The result is copied into a new AlgebraicMatrix.
   \return A new AlgebraicMatrix that is the transpose of the original.
   \sa AlgebraicMatrix<T>::transpose
 */
    // Allocate a matrix with ncols rows and nrows cols
    AlgebraicMatrix<T> temp(this->_ncols,this->_nrows);

    // Cycle over rows and colunms
    for(unsigned int i=0; i<(this->_nrows); i++)
        for(unsigned int j=0; j<(this->_ncols); j++)
            temp.at(j,i) = this->_data[i*(this->_ncols)+j];

    return temp;
}

/*! \cond */
/* Auxiliary functions for inverse and determinant computation, ignored in Doxygen Documentation */
template<class T> unsigned int AlgebraicMatrix<T>::pivot(unsigned int& k, const unsigned int ii, const AlgebraicMatrix<T>& A, std::vector<unsigned int>& P, std::vector<unsigned int>& R, std::vector<unsigned int>& C1, std::vector<unsigned int>& C2, T& det) {
    using std::abs;

    unsigned int im = 0;
    double t, m = 0;
    const unsigned int n = A.nrows();

    for (unsigned int i=0; i<n; i++){
        if (P[i]==0){
            for (unsigned int j=0; j<n; j++){
                if (P[j]==0){
                    t = abs( A.at(R[i],j) );
                    if (!(t<m)){
                        im = i;
                        k = j;
                        m = t;
                    }
                }
            }
        }
    }

    if (m<1.e-12){
        det = 0.0;
        return 123;}
    else{
        det = det*A.at( R[im], k ); // multiply pivot into determinant
        if (im!=k) det = -1.0*det;  // adjust sign of determinant
        std::swap(R[k], R[im]);
        P[k] = 1;                   // mark column as done

        // record original row/column of pivot
        C1[ii] = im;
        C2[ii] = k;

        return 0;}
}

/* Auxiliary functions for inverse and determinant computation */
template<class T> void AlgebraicMatrix<T>::eliminate( const unsigned int k, AlgebraicMatrix<T>& A, std::vector<unsigned int>& R) {
    unsigned int n = A.nrows();

//    // Better speed, less accuracy
//    for (unsigned int j=0; j<n; j++)
//    {
//        if (j!=k){
//            A.at(R[k],j) = A.at(R[k],j)*A.at(R[k],k);
//        }
//    }

    // Better accuracy, less speed
    for (unsigned int j=0; j<n; j++){
        if (j!=k)
            A.at(R[k],j) = A.at(R[k],j)/A.at(R[k],k);}

    A.at(R[k],k) = 1.0/A.at(R[k],k);

    for (unsigned int i=0; i<n; i++){
        if(i!=k){
            for (unsigned int j=0; j<n; j++){
                if (j!=k) A.at(R[i],j) = A.at(R[i],j) - A.at(R[i],k)*A.at(R[k],j);
            }
            A.at(R[i],k) = -1.0*A.at(R[i],k)*A.at(R[k],k);
        }
    }
}
/*! \endcond */

template<class T> AlgebraicMatrix<T> AlgebraicMatrix<T>::inv() const{
/*! Compute the inverse of an AlgebraicMatrix
   Algorithm based on the Gauss elimination with full pivot (from the Numerical Cookbook)
   The result is copied in a new AlgebraicMatrix.
   \return A new AlgebraicMatrix that is the inverse of the original.
   \sa AlgebraicMatrix<T>::inverse
 */
    // Check if matrix is square
    if (_nrows!=_ncols)
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::inv: Matrix must be square to compute the inverse.");

    const unsigned int n = _nrows;
    unsigned int k=0;
    T det = 1.0;
    std::vector<unsigned int> R(n);
    std::vector<unsigned int> P(n,0);
    std::vector<unsigned int> C1(n,0);
    std::vector<unsigned int> C2(n,0);
    AlgebraicMatrix<T> AI(n);
    AlgebraicMatrix<T> AA(*this);

    for (unsigned int i=0; i<n; i++) R[i] = i;

    for (unsigned int i=0; i<n; i++){
        unsigned int err = pivot(k, i, AA, P, R, C1, C2, det);
        if (err==0)
            eliminate(k, AA, R);
        else
            throw std::runtime_error("DACE::AlgebraicMatrix<T>::inv: Matrix inverse does not exist.");
    }

    for (unsigned int i=0; i<n; i++) P[i] = i;
    for (int i=n-1; i>=0; i--) std::swap(P[C1[i]], P[C2[i]]);
    for (unsigned int i=0; i<n; i++)
        for (unsigned int j=0; j<n; j++)
            AI.at(i,j) = AA.at(R[i], P[j]);

    return AI;
}

template<class T> T AlgebraicMatrix<T>::det() const{
/*! Compute the determinant of an AlgebraicMatrix
   \return A variable that is the determinant of the matrix.
   \sa AlgebraicMatrix<T>::det
 */
    // Check if matrix is square
    if (_nrows!=_ncols)
        throw std::runtime_error("DACE::AlgebraicMatrix<T>::inv: Matrix must be square to compute the inverse.");

    const unsigned int n = _nrows;
    unsigned int k=0;
    T det = 1.0;
    std::vector<unsigned int> R(n);
    std::vector<unsigned int> P(n,0);
    std::vector<unsigned int> C1(n,0);
    std::vector<unsigned int> C2(n,0);
    AlgebraicMatrix<T> AA(*this);

    for (unsigned int i=0; i<n; i++) R[i] = i;

    for (unsigned int i=0; i<n; i++){
        unsigned int err = pivot(k, i, AA, P, R, C1, C2, det);
        if (err==0)
            eliminate(k, AA, R);
        else
            return 0.0;
    }

    return det;
}

/***********************************************************************************
*     Coefficient access routines
************************************************************************************/
template<class T> AlgebraicMatrix<double> AlgebraicMatrix<T>::cons() const{
/*! Extract the constant part of a AlgebraicMatrix
   The result is copied in a new AlgebraicMatrix.
   \return A new AlgebraicMatrix that contains the constant part of the original.
   \sa AlgebraicMatrix<T>::cons
 */
    AlgebraicMatrix<double> temp(this->_nrows, this->_ncols);

    for(unsigned int i=0; i<(this->_nrows); i++)
        for(unsigned int j=0; j<(this->_ncols); j++)
            temp.at(i,j) = DACE::cons(this->_data[i*this->_ncols+j]);

    return temp;
}

/***********************************************************************************
*     Input/Output routines
************************************************************************************/
template<typename U> std::ostream& operator<< (std::ostream &out, const AlgebraicMatrix<U> &obj){
/*! Overload of std::operator<< in iostream.
   \param[in] out standard output stream.
   \param[in] obj AlgebraicMatrix<U> to be printed in the stream
   \return Standard output stream.
 */
    // Output each column at a time
    unsigned int nrows = obj.nrows();
    unsigned int ncols = obj.ncols();

    out << "[[[ " << nrows << "x" << ncols << " matrix" << std::endl;
    for(unsigned int i = 0; i<nrows;i++){
        for(unsigned int j = 0; j<ncols; j++){
            out << obj.at(i,j);
            if((j+1)==ncols)
                out<<std::endl;
            else
                out<<'\t';
        }
    }

  out << "]]]" << std::endl;

  return out;
}

template<typename U> std::istream& operator>> (std::istream &in, AlgebraicMatrix<U> &obj){
/*! Overload of std::operator>> in iostream.
    \param[in] in standard input stream.
    \param[in] obj AlgebraicMatrix<U> to be created from the stream
    \return Standard input stream.
 */
    // read the first line
    std::string init_line;
    std::getline(in, init_line);

    // initialize the size of the vector to be read
    unsigned int n_rows = 0;
    unsigned int n_cols = 0;

    // retrieve the size of the vector to be read
    if(in.good() ){
        // Find the number of rows
        std::size_t found = init_line.find_first_of("x");
        int j = 4;
        std::string size_str;
        while(j<(int)found){
            size_str = size_str + init_line[j];
            j++;}
        if (!(std::istringstream(size_str) >> n_rows)) n_rows = 0;

        // Find the number of columns (stop when matrix is met)
        j = (int)found+1;
        found = init_line.find_first_of("m",found);
        size_str.clear();
        while(j<(int)found){
            size_str = size_str + init_line[j];
            j++;}
        if (!(std::istringstream(size_str) >> n_cols)) n_cols = 0;

        // resize the object to meet the size of the vector to be read
        if ((obj.nrows() != n_rows)&&(obj.ncols() != n_cols))
            obj.resize(n_rows, n_cols);

        // fill the AlgebraicMatrix
        int i = 0;
        while(in.good() && (i<n_rows)) {
            // read value and allocate in the i-th component of obj
            for(int j=0; j<n_cols; j++)
                in >> obj.at(i,j);
            i++;
        }

        // check the next character
        if (in.peek() == '\n')    // the previous operator>> does not consume the \n character when an AlgebraicVector<T> (with T != DA) is considered
            in.ignore(); // ignore the next character

        // skip the line at the end of a AlgebraicMatrix (containing ]]])
        std::string skip_line;
        std::getline(in, skip_line);
    }
    return in;
}

/***********************************************************************************
*     Matrix operations
************************************************************************************/
template<class T> AlgebraicMatrix<T> transpose(const AlgebraicMatrix<T>& obj){
/*! Compute the transpose of an AlgebraicMatrix.
   \param[in] obj An AlgebraicMatrix
   \return An AlgebraicMatrix that is the transpose of the original AlgebraicMatrix.
   \sa AlgebraicMatrix<T>::transpose
 */
    return obj.transpose();
}

template<class T> T det(const AlgebraicMatrix<T>& obj){
/*! Compute the determinant of an AlgebraicMatrix.
   \param[in] obj An AlgebraicMatrix
   \return A variable that is the determinant of the AlgebraiMatrix.
   \sa AlgebraicMatrix<T>::det
 */
    return obj.det();
}

template<class T> AlgebraicMatrix<T> inv(const AlgebraicMatrix<T>& obj) {
/*! Compute the inverse of an AlgebraicMatrix.
   \param[in] obj An AlgebraicMatrix
   \return An AlgebraicMatrix that is the inverse of the original AlgebraicMatrix.
   \sa AlgebraicMatrix<T>::inverse
 */
    return obj.inv();
}

/***********************************************************************************
*     Coefficient access routines
************************************************************************************/
template<class T> AlgebraicMatrix<double> cons(const AlgebraicMatrix<T>& obj) {
/*! Extract the constant part of a AlgebraicMatrix
   The result is copied in a new AlgebraicMatrix.
   \param[in] obj An AlgebraicMatrix
   \return A new AlgebraicMatrix that contains the constant part of the original.
   \sa AlgebraicMatrix<T>::cons
 */
    return obj.cons();
}

}
#endif /* DINAMICA_DAMATRIX_T_H_ */
