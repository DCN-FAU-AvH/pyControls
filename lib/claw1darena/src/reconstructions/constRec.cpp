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

/*! \file constRec.cpp
 *  \brief Implementation of ConstRec
 */

#include "../config.h"
#include "../dof/dofvector.hpp"
#include "constRec.hpp"

//! \brief Constructor taking the grid in input
template <int M>
ConstRec<M>::ConstRec(Grid & grid):
  RecBase<M>(grid)
  {}

//! \brief Stores internally a pointer to the data being reconstructed
template <int M>
void ConstRec<M>::compute(const DofVec<M> * u/*! vector of conserved quantities associated to the grid*/, int, int){
  setUAvgAndResize(u);
}

//! Evaluation at any point: we just return the cell average
template <int M>
typename ConstRec<M>::t_u ConstRec<M>::eval(int k, tC){
  return (*_uData)[k];
}

template <int M>
const typename ConstRec<M>::t_u * ConstRec<M>::getCoeffRec(int k) const
{return 0;}

template <int M>
void ConstRec<M>::overrideCoeffRec(int k, const t_u* coeff)
{}

//! \brief Number of ghost cells (0)
template <int M>
int ConstRec<M>::needsGhosts() const{
  /*! We don't need ghost cells for the reconstruction, thus return 0 */
  return 0;
}

//! \brief Order of accuracy (1)
template <int M>
int ConstRec<M>::getOrder() const{
  return 1;
}

template
class ConstRec<1>;

template
class ConstRec<2>;

template
class ConstRec<3>;

/** @} */
