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

/** @addtogroup timeintegration
 *
 * @{
 */

/*! \file timeEuler.tpp
 *  \brief Implementation of Euler
 */

#include "assert.h"

#include "timeEuler.hpp"
#include "../grid/gridutils.hpp"

template<int M>
Euler<M>::Euler(Grid & grid, const semidiscreteRHS<M> & doRHS):
//RecBase<M> & rec, FluxBase<M> & numFlux, BCHandler<M>& bcHandler, numSource<M>& source):
  //timeIntegration<M>(grid, rec, numFlux, bcHandler, source)
  timeIntegration<M>(grid),
  _doRHS(doRHS)
  {}

/*! \brief Explicit Euler timestepping scheme
 *
 * U0 is assumed to contain data associated to the grid passed
 * to the constructor.
 */
template<int M>
tC Euler<M>::advance(t_U &U0, tC t0, t_U &U1){

  /*! - U1 is resized to the size of U0 */
  U1.resize(U0.size());

  ///*! - the fluxes, sources and thus the stage K is computed. */
  tD vMax(0.);
  _doRHS.computeRHS(U0, t0, K, vMax);

  /*! - The timestep dt is the one stored in _dt if _useCFL is false, and otherwise is computed as _cfl*dx/vMax, where vMax is the largest wave speed encountered during the flux computation.
   */
  const tC dx = timeIntegration<M>::_grid.getDx();
  if (timeIntegration<M>::_useCFL){
    timeIntegration<M>::_dt = timeIntegration<M>::_cfl*dx/vMax;
    if (timeIntegration<M>::_dt>timeIntegration<M>::_dtMax)
      timeIntegration<M>::_dt=timeIntegration<M>::_dtMax;
  }
  const tD lambda(timeIntegration<M>::_dt/dx); //conversion of dt/dx to type tD

  /*! - Finally, U1 is computed in the physical cells (begin/endPhysical
   *  Grid methods)
   *
   */

  for (int i=timeIntegration<M>::_grid.beginPhysical(); i<timeIntegration<M>::_grid.endPhysical(); ++i){
    U1[i] = U0[i] + lambda*K[i];
    //TODO: cambiare il segno di K per fare K=-deltaF+sorgente e mettere +lambda*K qui
  }

  /*!
   * Upon exit, U1 will contain the updated solution (only physical
   * cells) and the timestep employed is returned
   */
  return (timeIntegration<M>::_dt);
}

//! \brief Number of ghost cells (1)
template <int M>
int Euler<M>::needsGhosts() const{
  return _doRHS.needsGhosts();
}

//! \brief Order of accuracy (1)
template <int M>
int Euler<M>::getOrder() const{
  return 1;
}

/** @} */
