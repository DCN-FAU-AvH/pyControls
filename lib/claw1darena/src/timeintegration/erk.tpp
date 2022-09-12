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

/*! \file erk.hpp
 *  \brief Declaration of ExplicitRungeKutta
 */

//#include "assert.h"

#include "erk.hpp"
//#include "../grid/gridutils.hpp"

//! \brief Constructor from tableaux, grid, reconstruction, flux and bcHandler
template<int M>
explicitRungeKutta<M>::explicitRungeKutta(ExplicitButcherTableaux & tableaux, Grid & grid, const semidiscreteRHS<M> & doRHS):
  timeIntegration<M>(grid),
  _doRHS(doRHS),
  _tableaux(tableaux)
  {
    _nStages = _tableaux.getNStages();
    _K = new t_U[_nStages];
  }

//! \brief Destructor: free space allocated for the stages
template<int M>
explicitRungeKutta<M>::~explicitRungeKutta(){
  delete[] _K;
}

//! \brief We need 1 ghost to compute the flux across the physical domain boundary
template<int M>
int explicitRungeKutta<M>::needsGhosts() const{
  return _doRHS.needsGhosts();
}

/*! \brief Explicit Runge-Kutta timestepping scheme
 *
 * U0 is assumed to contain data associated to the grid passed
 * to the constructor.
 */
template<int M>
tC explicitRungeKutta<M>::advance(t_U &U0, tC t0, t_U &U1){

  /*! - U1 is resized to the size of U0 */
  U1.resize(U0.size());
  const tC dx = timeIntegration<M>::_grid.getDx();
  tD lambda(0.);

  /*! - Loop on stages: */
  for (int k=0; k<_nStages; ++k){

    /*!   * compute the stage value usign the A_{ij} coefficients of the Runge-Kutta*/
    for (int i=timeIntegration<M>::_grid.beginPhysical(); i<timeIntegration<M>::_grid.endPhysical(); ++i){
      U1[i] = U0[i];
      // In the first stage, lambda=0, but we won't enter the following loop
      for (int j=0; j<k; ++j)
        U1[i] += lambda * _tableaux.getA(k,j) * _K[j][i];
    }

    /*!   * set the boundary conditions */
    // In the first stage _dt is not set yet, but getC(0)=0
    const tC t1=t0 + _tableaux.getC(k)*timeIntegration<M>::_dt;

    /*! * Compute the r.h.s. of the semidiscrete equation on the current stage value */
    tD vMax(0.);
    _doRHS.computeRHS(U1, t0, _K[k], vMax);

    /*!  * (only in the first stage) The timestep dt is the one stored in _dt if _useCFL is false, and otherwise is computed as _cfl*dx/vMax, where vMax is the largest wave speed encountered during the flux computation.
     */
    if ((timeIntegration<M>::_useCFL) && (k==0)){
      timeIntegration<M>::_dt = timeIntegration<M>::_cfl*dx/vMax;
      if (timeIntegration<M>::_dt>timeIntegration<M>::_dtMax)
        timeIntegration<M>::_dt=timeIntegration<M>::_dtMax;      
    }
    //AS: bugfix, lambda has to be set also for a fixed dt, as dt may change in
    // each timestep
    lambda = tD(timeIntegration<M>::_dt/dx); //conversion of dt/dx to type tD
  }

  /*! - The final solution is computed in the physical cells (begin/endPhysical
   *  Grid methods) using the b_j coefficients of the Runge-Kutta
   */
  for (int i=timeIntegration<M>::_grid.beginPhysical(); i<timeIntegration<M>::_grid.endPhysical(); ++i){
    U1[i] = U0[i];
    for (int j=0; j<_nStages; ++j)
      U1[i] += lambda * _tableaux.getB(j) * _K[j][i];
  }

  /*!
   * Upon exit, U1 will contain the updated solution (only physical
   * cells) and the timestep employed is returned
   */
  return (timeIntegration<M>::_dt);
}

//! \brief Order of accuracy: returns the order reported by the ExplicitButcherTableaux
template<int M>
int explicitRungeKutta<M>::getOrder() const{
  return _tableaux.getOrder();
}

/** @} */
