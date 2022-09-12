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

/*! \file bcHandler.hpp
 *  \brief Declaration of BCHandler
 */

#ifndef _BCHANDLERBASE_HH
#define _BCHANDLERBASE_HH

#include "../grid/grid.hpp"
#include "../dof/dofvector.hpp"

/*! \brief Base class for handling boundary conditions
 *
 * Each Conservation law implementation should also provide a derived
 * class CLawNameBCHandler that handles the boundary conditions that
 * are significant for that specific conservation law.
 *
 * This is just a pure virtual  base class.
 */
template <int M>
class BCHandler{
public:

  virtual ~BCHandler() {} //!< virtual destructor

  /*! \brief Sets the values of the left ghosts in U
   *
   * In input, U is assumed to be a vector of conserved quantities
   * associated to the grid.
   * In output, the elements of U corresponding to the left ghosts are
   * set according to the specified boundary condition and the data in
   * the physical cells
   */
  virtual void setLeftGhosts  (Grid &grid /*! grid*/, Grid::BCType bc /*! boundary condition type*/, tC t /*!current time*/, DofVec<M> &U /*! vector of conserved variables on the grid*/)=0;

  /*! \brief Sets the values of the left ghosts in U
   *
   * In input, U is assumed to be a vector of conserved quantities
   * associated to the grid.
   * In output, the elements of U corresponding to the right ghosts are
   * set according to the specified boundary condition and the data in
   * the physical cells
   */
  virtual void setRightGhosts (Grid &grid /*! grid*/, Grid::BCType bc /*! boundary condition type*/, tC t /*!current time*/, DofVec<M> &U /*! vector of conserved variables on the grid*/)=0;
};

//! \brief Base functor class for time dependent boundary data
template <int M>
class bcFuncType{
public:
  typedef DoF<M,tD> t_u;
  //! \brief operator() that returns the fuction value
  virtual t_u operator() (tC t);
};

#endif

/** @} */
