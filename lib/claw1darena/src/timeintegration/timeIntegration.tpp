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

/** @defgroup timeintegration Time Integration
 *
 * @{
 */

/*! \file timeIntegration.tpp
 *  \brief Implementation of class naifRHS
 */

#include "timeIntegration.hpp"

template <int M>
void naifRHS<M>::computeRHS(t_U & U, tD time, t_U &rhs, tD &vMax) const{

    /*! - First, the boundary conditions on U are set */
  _bcHandler.setLeftGhosts (_grid , _grid.getLeftBC() , time, U);
  _bcHandler.setRightGhosts(_grid , _grid.getRightBC(), time, U);

  /*! - rhs is resized to the size of U */
  rhs.resize(U.size());
  rhs.zero();

  /*! - the reconstruction of the data in U is computed */
  _rec.compute(&U);
  vMax=0;

  //! - compute the flux contribution of each face (using the Grid methods begin/endFaces, left/rightOfFace and isPhysical) and stack them in rhs
  for (int f = _grid.beginFaces(); f<_grid.endFaces(); ++f){
    int leftId = _grid.leftOfFace(f);
    t_u uL = _rec.evalRight(leftId);
    int rightId = _grid.rightOfFace(f);
    t_u uR = _rec.evalLeft(rightId);
    tD vmaxL, vmaxR;
    t_u F = _numFlux.getF(uL,uR,0,0,vmaxL,vmaxR);
    vMax=(vmaxL>vMax? vmaxL: vMax);
    vMax=(vmaxR>vMax? vmaxR: vMax);
    if (_grid.isPhysical(leftId))
      rhs[leftId] -= F;
    if (_grid.isPhysical(rightId))
      rhs[rightId] += F;
  }

  //! - add the source in the physical cells (begin/endPhysical Grid methods)
  for (int i=_grid.beginPhysical(); i<_grid.endPhysical(); ++i){
    rhs[i] += tD(_grid.getDx()) * _source.getSource(i,time);
  }
}

template <int M>
int naifRHS<M>::needsGhosts() const{
  return _rec.needsGhosts() + 1;
}

template <int M>
int naifRHS<M>::getOrder() const{
  return _rec.getOrder();
}


/** @} */
