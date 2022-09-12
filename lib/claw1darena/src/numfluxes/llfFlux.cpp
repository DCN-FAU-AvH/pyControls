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

/*! \file llfFlux.cpp
 *  \brief Implementation of llfFlux
 */

#include "llfFlux.hpp"
#include "../claws/claw_linear.hpp"
#include "../claws/claw_burgers.hpp"
#include "../claws/claw_gas.hpp"
#include "../claws/claw_psystem.hpp"
#include "../claws/claw_swe.hpp"

/*! \brief Constructor
 *
 * Pass in input an object of class CLaw.
 *
 * This will be the llfFlux computational workhorse: it will be called
 * whenever an exact flux or a spectral radius of a flux Jacobian is
 * needed
 */
template <class CLaw>
llfFlux<CLaw>::llfFlux(CLawBase<CLaw,CLaw::m> & cLaw):
  _cLaw(cLaw)
  {}

//! \brief Computes the flux
template <class CLaw>
typename llfFlux<CLaw>::t_u llfFlux<CLaw>::getF(const t_u& uL/*!left state*/, const t_u& uR/*!right state*/, tC x /*! position of the face*/, tC t /*! time*/, tD & lL /*! set to &alpha; in output*/, tD& lR/*! set to &alpha; in output*/)
{
  t_u F;
  F = _cLaw.getFlMax(uL,x,t,lL);
  F += _cLaw.getFlMax(uR,x,t,lR);
  if (lL<lR)
    lL=lR;
  else
    lR=lL;
  F -= lL * (uR-uL);
  F *= 0.5;
  return F;
}

// Generate code for the specialization for CLawLinear
template
class llfFlux<CLawLinear>;

// Generate code for the specialization for CLawBurgers
template
class llfFlux<CLawBurgers>;

// Generate code for the specialization for CLawGas
template
class llfFlux<CLawGas>;

// Generate code for the specialization for CLawGas
template
class llfFlux<CLawPsys>;


// Generate code for the specialization for CLawSWE

//! \brief Computes the flux
template<>
CLawSWE::t_u llfFlux<CLawSWE>::getF(const CLawSWE::t_u& uL/*!left state*/, const CLawSWE::t_u& uR/*!right state*/, tC x /*! position of the face*/, tC t /*! time*/, tD & lL /*! set to &alpha; in output*/, tD& lR/*! set to &alpha; in output*/)
{
  //usual code for LaxFriedrichs
  t_u F;
  F = _cLaw.getFlMax(uL,x,t,lL);
  F += _cLaw.getFlMax(uR,x,t,lR);
  if (lL<lR)
    lL=lR;
  else
    lR=lL;
  F -= lL * (uR-uL);
  F *= 0.5;
  //but no flux for Z variable
  F[CLawSWE::Z]=0.;
  return F;
}

template
class llfFlux<CLawSWE>;

/** @} */
