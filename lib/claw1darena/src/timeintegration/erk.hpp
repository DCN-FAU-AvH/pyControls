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

#ifndef ERK_HH
#define ERK_HH

#include "timeIntegration.hpp"
#include "butcher.hpp"

/*! \brief Explicit Runge-Kutta time integration
 *
 * This class implements a time-integration procedure based on an
 * explicit Runge-Kutta scheme. Any explicit Runge-Kutta tableaux
 * can be specified via the constructor of the class.
 *
 * The template parameter is for the number of conservation laws in the
 * system of equations.
 *
 * The r.h.s. of the equation is computed with the semidiscreteRHS object passed to the constructor
 */
template<int M>
class explicitRungeKutta: public timeIntegration<M>{
public:
  typedef typename timeIntegration<M>::t_U t_U; //!< type for a vector of conserved quantities on the grid

  explicitRungeKutta(ExplicitButcherTableaux & tableaux,Grid & grid, const semidiscreteRHS<M> & doRHS);

  ~explicitRungeKutta();

  virtual int needsGhosts() const;

  virtual tC advance(t_U &u0, tC t0, t_U &u1);

  virtual int getOrder() const;
private:
  ExplicitButcherTableaux & _tableaux;
  int _nStages;
  DofVec<M> *_K; //!< Internal storage for the stages
  const semidiscreteRHS<M> &_doRHS; //!< Reference to the object computing the r.h.s. of the system
};

#endif

/** @} */
