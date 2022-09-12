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

/*! \file timeEuler.hpp
 *  \brief Declaration of Euler
 */

#ifndef TIMEEULER_HH
#define TIMEEULER_HH

#include "timeIntegration.hpp"

/*! \brief Explicit Euler time integration
 * *
 * The template parameter is for the number of conservation laws in the
 * system of equations.
 *
 * The r.h.s. of the equation is computed with the semidiscreteRHS object passed to the constructor
 */
template<int M>
class Euler: public timeIntegration<M>{
public:
  //! \brief Constructor: sets grid, reconstruction, numerical flux and boundary conditions
  Euler(Grid & grid, const semidiscreteRHS<M> & doRHS);
  // RecBase<M> & rec, FluxBase<M> & numFlux, BCHandler<M>& bcHandler, numSource<M>& source);

  //typedef typename timeIntegration<M>::t_u t_u; //!< type for a set of conserved quantities
  typedef typename timeIntegration<M>::t_U t_U; //!< type for a vector of conserved quantities on the grid

  virtual tC advance(t_U &u0, tC t0, t_U &u1);

  virtual int needsGhosts() const ;

  virtual int getOrder() const;

private:
  DofVec<M> K; //!< Internal storage for the stage
  const semidiscreteRHS<M> &_doRHS; //!< Reference to the object computing the r.h.s. of the system
};

#endif

/** @} */
