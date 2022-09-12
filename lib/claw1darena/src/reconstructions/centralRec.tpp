/*
  This file is part of claw1dArena, a software library to discretize
  one-dimensional conservation laws.

  Copyright (C) 2017 Matteo Semplice

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

  For details on the licence, see the COPYING file at the root of the tree.
*/

/** @addtogroup reconstructions
 *
 * @{
 */

/*! \file centralRec.tpp
 *  \brief Implementation of CentralRec
 */

#ifndef CENTRALREC_TPP
#define CENTRALREC_TPP

#include "centralRec.hpp"
#include "udivdiff.hpp"
#include "polyTools.hpp"

//! \brief Constructor from grid and accuracy
template <int M>
CentralRec<M>::CentralRec(Grid & grid, int accuracy):
  RecBase<M>(grid),
  _udd()
  {setAccuracy(accuracy);}

/*! \brief Sets the accuracy of the reconstruction
 *
 * Note that, for the reconstruction to be central, the accuracy has to
 * be an odd integer. This in enforced here by increasing the accuracy
 * by 1 if an even integer is passed.
 *
 * This method also resizes the internal storage variables.
 */
template <int M>
void CentralRec<M>::setAccuracy(int accuracy){
  _accuracy = 2*(accuracy/2)+1; //makes sure its odd...
}

//! \brief Returns the current accuracy of the reconstruction
template <int M>
int  CentralRec<M>::getAccuracy(){
  return _accuracy;
}

//! \brief Computes the reconstruction polynomials.
template <int M>
void CentralRec<M>::compute(const DofVec<M> * u, int cellStart, int cellEnd){
  setUAvgAndResize(u);
  _udd.compute(*u,_accuracy-1);
  int r = _accuracy/2; //shift for central polynomial
  for (int k=cellStart; k<cellEnd; ++k){
    _udd.InterpPoly(k, _accuracy-1, r, _poly[k].data() );
  }
}

/*! \brief evaluates the reconstruction at position pos
 *
 * pos is the relative position within the cell, thus it should be
 * between -0.5 and 0.5.
 */
template <int M>
typename CentralRec<M>::t_u CentralRec<M>::eval(int k, tC pos){
  return evalPoly(_accuracy-1, pos, _poly[k].data() );
}

template <int M>
const typename CentralRec<M>::t_u * CentralRec<M>::getCoeffRec(int k) const{
  return _poly[k].data();
}

template <int M>
void CentralRec<M>::overrideCoeffRec(int k, const t_u* coeff){
  for (int d=0; d<_accuracy+1; ++d)
    _poly[k][d] = coeff[d];
}

//! \brief Returns the number of ghost cells used by this reconstruction
template <int M>
int CentralRec<M>::needsGhosts() const{
  return _accuracy/2;
}

//! \brief Returns the order of accuracy
template <int M>
int CentralRec<M>::getOrder() const{
  return _accuracy;
}

//! \brief Resize storage for reconstruction polynomial
template <int M>
void CentralRec<M>::setUAvgAndResize(const DofVec<M> * u){
  if ( (_poly.size() < u->size() ) || (_poly[0].size() < _accuracy) ) {
    _poly.resize(u->size());
    for (int k=0; k<u->size(); ++k)
      _poly[k].resize(_accuracy+1);
  }
}

#endif

/** @} */
