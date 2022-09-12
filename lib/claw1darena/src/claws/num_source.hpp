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

/** @addtogroup claws
 *
 * @{
 */

/*! \file num_source.hpp
 *  \brief Declaration of numSource
 */

#ifndef _NUMSOURCE_HH
#define _NUMSOURCE_HH

#include "../grid/grid.hpp"

/*! \brief Base class for numerical computation of source terms
 *
 * The implementation of gas dynamics should also provide a derived
 * class for each source term that handles the quadrature of the source
 * term for that each specific case.
 *
 * This is just a pure virtual base class.
 */
template <int M>
class numSource{
public:

  //! \brief Default constructor
  numSource(Grid & grid):
    _grid(grid)
    {}

  virtual ~numSource() {} //!< virtual destructor

  typedef DoF<M,tD> t_u; //!< type for a set of conserved quantities

  /*! \brief Computes the source cell average
   *
   * Purely virtual since it can only be implemented in the claw*.tpp files
   */
  virtual t_u getSource(int k /*! cell index*/, tC t /*! time*/) const =0;

protected:
    Grid & _grid; //!< reference to the grid
};

template <int M>
class ZeroSource : public numSource<M>{
public:
  typedef DoF<M,tD> t_u; //!< type for a set of conserved quantities

  ZeroSource(Grid & grid):
    numSource<M>(grid)
  {}

  //! \brief returns the identically zero source term for a system of M equations
  virtual t_u getSource(int k, tC t) const {
    return t_u(0.0);
  }
};


#endif

/** @} */
