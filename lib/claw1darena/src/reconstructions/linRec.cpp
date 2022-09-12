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

/*! \file linRec.cpp
 *  \brief Implementation of LinRec
 */

#include "../config.h"
#include "../dof/dof.hpp"
#include "../dof/dofvector.hpp"
#include "linRec.hpp"
#include <cmath>

template <int M>
LinRec<M>::LinRec(Grid & grid):
  RecBase<M>(grid),
  _slopeData()
  {}

/*! \brief Computes the slope-limited linear reconstruction
 *
 * Unlimited slopes are first computed and then the minmod limiter is
 * applied to compute the limited slopes, that are stored in _slopeData.
 */
template <int M>
void LinRec<M>::compute(const DofVec<M> * u, int cellStart, int cellEnd){
  setUAvgAndResize(u);
  _uSlope.resize(u->size());

  for(int i=cellStart-1; i<cellEnd; ++i)
    _uSlope[i] = (*u)[i+1] -(*u)[i];

  for (int i=cellStart; i<cellEnd; ++i){
    _slopeData[i]=t_u(0.);

    const t_u & lSlope = _uSlope[i-1];
    const t_u & rSlope = _uSlope[i];

    for (int k=0; k<M; ++k){
      if ( lSlope[k] * rSlope[k] <0.)
        continue;
      if (std::abs(lSlope[k]) < std::abs(rSlope[k]) )
        _slopeData[i][k] = lSlope[k];
      else
        _slopeData[i][k] = rSlope[k];
    }
  }
}

template <int M>
typename LinRec<M>::t_u LinRec<M>::eval(int k, tC relx){
  return ((*_uData)[k] + tD(relx)*(_slopeData)[k] );
}

template <int M>
const typename LinRec<M>::t_u * LinRec<M>::getCoeffRec(int k) const{
  return &(_slopeData[k]);
}

template <int M>
void LinRec<M>::overrideCoeffRec(int k, const t_u* coeff){
  _slopeData[k] = coeff[0];
}

//! \brief Number of ghost cells (1)
template <int M>
int LinRec<M>::needsGhosts() const{
  /*! We need 1 ghost cell per side to compute the reconstruction, thus return 1 */
  return 1;
}

//! \brief Sets pointer to data being reconstructed and resize storage for limited slopes
template <int M>
void LinRec<M>::setUAvgAndResize(const DofVec<M> * u){
  _uData=u; 
  _slopeData.resize(u->size());
}

//! \brief Order of accuracy (2)
template <int M>
int LinRec<M>::getOrder() const{
  return 2;
}

template
class LinRec<1>;

template
class LinRec<3>;

/** @} */
