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

/** @addtogroup grid
 *
 * @{
 */

/*! \file gridutils.cpp
 *  \brief Definitions of helper functions to set ghost values
 */

#include "gridutils.hpp"

/*! \brief set left ghosts for periodic boundary conditions */
template <int M>
void setLeftPeriodicGhosts(Grid &G , DofVec<M> &U){
  int k=G.endPhysical()-1; //last physical cell
  for (int i=G.rbeginLeftGhosts(); i>G.rendLeftGhosts() ; i--){
    U[i] = U[k];
    k--;
  }
}

/*! \brief set right ghosts for periodic boundary conditions */
template <int M>
void setRightPeriodicGhosts(Grid &G , DofVec<M> &U){
  int k=G.beginPhysical();
  for (int i=G.beginRightGhosts(); i<G.endRightGhosts() ; i++){
    U[i] = U[k];
    k++;
  }
}

/*! \brief set left ghosts for free-flow boundary conditions */
template <int M>
void setLeftFreeFlowGhosts(Grid &G , DofVec<M> &U){
  int k=G.beginPhysical(); //first physical cell
  for (int i=G.rbeginLeftGhosts(); i>G.rendLeftGhosts() ; i--)
    U[i] = U[k];
}

/*! \brief set right ghosts for free-flow boundary conditions */
template <int M>
void setRightFreeFlowGhosts(Grid &G , DofVec<M> &U){
  int k=G.endPhysical()-1; //last physical cell
  for (int i=G.beginRightGhosts(); i<G.endRightGhosts() ; i++)
    U[i] = U[k];
}

/*! \brief set left ghosts for Dirichlet boundary conditions */
template <int M>
void setLeftDirichletGhosts(Grid &G , DofVec<M> &U, DoF<M,tD> value){
  int k=G.beginPhysical(); //first physical cell
  int d=1;
  for (int i=G.rbeginLeftGhosts(); i>G.rendLeftGhosts() ; i--, d++)
    U[i] = tD(2*d)*value - tD(2*d-1)*U[k];
}

/*! \brief set right ghosts for Dirichlet boundary conditions */
template <int M>
void setRightDirichletGhosts(Grid &G , DofVec<M> &U, DoF<M,tD> value){
  int k=G.endPhysical()-1; //last physical cell
  int d=1;
  for (int i=G.beginRightGhosts(); i<G.endRightGhosts() ; i++, d++)
    U[i] = tD(2*d)*value - tD(2*d-1)*U[k];
}

template
void setLeftPeriodicGhosts<1>(Grid &G , DofVec< 1 > &U);
template
void setRightPeriodicGhosts<1>(Grid &G , DofVec< 1 > &U);

template
void setLeftFreeFlowGhosts<1>(Grid &G , DofVec< 1 > &U);
template
void setRightFreeFlowGhosts<1>(Grid &G , DofVec< 1 > &U);

template
void setLeftDirichletGhosts<1>(Grid &G , DofVec< 1 > &U, DoF<1,tD> value);
template
void setRightDirichletGhosts<1>(Grid &G , DofVec< 1 > &U, DoF<1,tD> value);

template
void setLeftPeriodicGhosts<3>(Grid &G , DofVec< 3 > &U);
template
void setRightPeriodicGhosts<3>(Grid &G , DofVec< 3 > &U);

template
void setLeftFreeFlowGhosts(Grid &G , DofVec< 3 > &U);
template
void setRightFreeFlowGhosts(Grid &G , DofVec< 3 > &U);

template
void setLeftDirichletGhosts(Grid &G , DofVec< 3 > &U, DoF<3,tD> value);
template
void setRightDirichletGhosts(Grid &G , DofVec< 3 > &U, DoF<3,tD> value);

/** @} */
