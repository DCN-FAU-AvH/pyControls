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

/** @addtogroup numfluxes
 *
 * @{
 */

/*! \file llfFlux.hpp
 *  \brief Declaration of llfFlux
 */

#ifndef LLFFLUX_HH
#define LLFFLUX_HH

#include "../config.h"
#include "../dof/dof.hpp"
#include "../claws/clawbase.hpp"
#include "fluxBase.hpp"

/*! \brief Class for the Local Lax Friedrichs (aka Rusanov) numerical flux
 *
 * The class is templated on the CLaw (which should be a derived class
 * of CLawBase)
 *
 * The numerical flux is given by
 *
 * F(uL,uR) = 0.5(f(uL)+f(uR)) - 0.5&alpha;(uR-uL)
 *
 * with &alpha; set to the largest of the spectral radii of f'(uL) and f'(uR).
 *
 * An object of type CLaw is passed to the constructor and a reference is
 * stored internally, so that it can be used to compute the eaxct flux
 * on a state and the spectral radii of the flux Jacobian.
 */
template <class CLaw>
class llfFlux: public FluxBase<CLaw::m> {
public:
  typedef typename CLaw::t_u t_u;          //!< type for a set of conserved quantities
  typedef CLawBase<CLaw,CLaw::m> cLaw_t; //!< type for the CLaw object

  llfFlux(cLaw_t & cLaw);

  virtual t_u getF(const t_u& uL, const t_u& uR, tC x, tC t, tD & lL, tD& lR);
private:
  cLaw_t & _cLaw; //!< reference to an object of class CLaw
};

#endif

/** @} */
