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

/** @defgroup numfluxes Numerical fluxes
 *
 * @{
 */

/*! \file fluxBase.hpp
 *  \brief Declaration of FluxBase
 */

#ifndef FLUXBASE_HH
#define FLUXBASE_HH

/*! \brief Base class for numerical fluxes
 *
 * Pure virtual class, templated on the number of conservative variables.
 *
 * In order to define a new numerical flux, derive a class from this one
 * and implement the getF method.
 */
template <int M>
class FluxBase{
public:
  typedef DoF<M, tD> t_u; //!< type for a set of conserved quantities

  //! \brief Default constructor
  FluxBase()
    {}

  virtual ~FluxBase() {} //!< virtual destructor

  /*! \brief Computes the numerical flux
   *
   * Returns the numerical flux at a given location and time, for the
   * left and right states passed in input.
   * Each implementation is also expected to return (an estimate for)
   * the maximal speed of the left-going and of the right-going waves.
   */
  virtual t_u getF(const t_u& uL/*!left state*/, const t_u& uR/*!right state*/, tC x/*! location*/, tC t/*! time*/, tD & lL/*! maximal speed of left-going waves (output)*/, tD& lR/*! maximal speed of left-going waves (output)*/) =0;
};

#endif

/** @} */
