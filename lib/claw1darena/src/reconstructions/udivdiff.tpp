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

/*! \file udivdiff.hpp
 *  \brief Declaration of UDivDiff
 */

#ifndef UDIVDIFF_TPP
#define UDIVDIFF_TPP

#include "udivdiff.hpp"

//! \brief Default constructor
template <int M>
UDivDiff<M>::UDivDiff()
  {}

/*! \brief Compute the (un)divided differences up to level maxD
 *
 * The table is stored internally in _udd and a pointer to the cell
 * averages u is stored internally and is used later by InterPoly: the
 * data should not change between the call to this function and the
 * call to InterPoly!
 */
template <int M>
void UDivDiff<M>::compute(DofVec<M> const & u/*! the vector of cell averages of conserved variables*/,int maxD){
  _udd.resize(u.size());
  _u=&u;
  for(int i=0; i<u.size()-1; ++i){
    _udd[i].resize(maxD);
    _udd[i][0] = tD(0.5) * (u[i+1] - u[i]);
  }
  for(int j=1; j<maxD; ++j){
    for(int i=0; i<u.size()-1-j; ++i){
      _udd[i][j] = (_udd[i+1][j-1] - _udd[i][j-1]) / (j+tD(2.));
    }
}

}

/*! Returns the (un)divided difference of order deg for the given cell
 *
 * Warning: no range-checking is made!
 * In particular, the highest degree (un)divided differences for the last
 * cells in the grid do not exist!
 */
template <int M>
typename UDivDiff<M>::t_u UDivDiff<M>::getUDD(int cell, int deg){
  if(deg>0)
    return _udd[cell][deg-1];
  else
    return (*_u)[cell];
}

/*! Returns the coefficients of the interpolating polynomial
 *
 * for a polynomial of degree deg, shift chooses the stencil among the
 * deg+1 possible stencils that include cell. It should be between 0
 * (left-most possible stencil) and deg (right-most possible stencil).
 *
 * In input, P should have space allocated for deg+1 coefficients.
 *
 * In output, the coefficients of the polynomial are stored in P,
 * starting from the least significant one.
 */
template <int M>
void UDivDiff<M>::InterpPoly(int cell, int deg, int shift,t_u* P){
  for (int j=0;j<deg+1;++j)
    P[j]=t_u(0.0);

  int row=cell-deg+shift;
  tD xf=shift-0.5;
  t_u tn=t_u(0);

  // primitive = tn + P[0]*x + P[1]*x^2 + ...
  for (int d=deg;d>=0;--d){
    tn+=getUDD(row,d);
    // multiply by x-xf
    for (int i=deg;i>0;--i)
      P[i]=(-xf)*P[i]+P[i-1];
    P[0]=(-xf)*P[0]+tn;
    tn*=-xf;
    xf-=1;
  }
  // differentiate P
  for (int k=0;k<deg+1;++k)
    P[k]*=(k+1);
}

#endif

/** @} */
