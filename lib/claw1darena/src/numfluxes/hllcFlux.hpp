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

/*! \file hllcFlux.hpp
 *  \brief Declaration of hllcFlux
 */

#ifndef HLLC_HH
#define HLLC_HH

#include "../config.h"
#include "../dof/dof.hpp"
#include "../claws/clawbase.hpp"
#include "../claws/claw_linear.hpp"
#include "../claws/claw_burgers.hpp"
#include "../claws/claw_gas.hpp"
#include "fluxBase.hpp"

/*! \brief Base class for the HLLC numerical flux
 *
 * HLLC is a class of approximate Riemann solvers that use 2
 * intermediate states (i.e. 3 waves) and for which the middle wave
 * is a contact. See Chap 10 of Toro's book for details.
 *
 * This a pure virtual class, templated on CLaw (which should be derived
 * from CLawBase). The derived class will provide the actual
 * implementations.
 */
template <class CLaw>
class hllcFlux: public FluxBase<CLaw::m> {
public:
  typedef DoF<CLaw::m, tD> t_u; //!< type for a set of conserved quantities
  typedef CLawBase<CLaw,CLaw::m> cLaw_t;

  /*! \brief Constructor
   *
   * Pass in an object of class CLaw: it will be used when fluxes or
   * spectral radii of flux Jacobians are needed.
   */
  hllcFlux(cLaw_t & cLaw);

  /*! \brief Returns the numerical flux
   *
   * In input, the states and the spatial and time location are passed.
   *
   * In output, the flux is returned and also the speed of the first and
   * third vawes are returned.
   */
  virtual t_u getF(const t_u& uL, const t_u& uR, tC x, tC t, tD & lL, tD& lR)=0;
private:
  cLaw_t & _cLaw;
};

template <>
class hllcFlux<CLawLinear>: public FluxBase<CLawLinear::m> {
public:
  typedef DoF<CLawLinear::m, tD> t_u;
  typedef CLawLinear cLaw_t;

  hllcFlux(CLawLinear & cLaw);

  virtual t_u getF(const t_u& uL, const t_u& uR, tC x, tC t, tD & lL, tD& lR);
private:
  cLaw_t & _cLaw;
};

template <>
class hllcFlux<CLawBurgers>: public FluxBase<CLawBurgers::m> {
public:
  typedef DoF<CLawBurgers::m, tD> t_u;
  typedef CLawBurgers cLaw_t;

  hllcFlux(CLawBurgers & cLaw);

  virtual t_u getF(const t_u& uL, const t_u& uR, tC x, tC t, tD & lL, tD& lR);
private:
  cLaw_t & _cLaw;
};

template <>
class hllcFlux<CLawGas>: public FluxBase<CLawGas::m> {
public:
  typedef DoF<CLawGas::m, tD> t_u;
  typedef CLawGas cLaw_t;

  hllcFlux(CLawGas & cLaw);

  virtual t_u getF(const t_u& uL, const t_u& uR, tC x, tC t, tD & lL, tD& lR);
private:
  cLaw_t & _cLaw;
};

#endif

/** @} */
