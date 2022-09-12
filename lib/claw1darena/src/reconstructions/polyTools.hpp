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

/** @addtogroup reconstructions
 *
 * @{
 */

/*! \file polyTools.hpp
 *  \brief Declaration of polyTools fucntions
 */

#ifndef POLYTOOLS_HH
#define POLYTOOLS_HH

#include "../grid/grid.hpp"
#include "../dof/dof.hpp"
#include "../dof/dofvector.hpp"
#include <vector>
#include <array>

/*! \brief Tools for polynomials
 *
 * Functions for computing indicators, evaluating a polynomials etc 
 * 
 * All the functions are templated on the number of conserved
 * quantities M.
 */

template<int M>
DoF<M,tD> computeIndicator(int deg, DoF<M,tD> *P);

template<int M>
DoF<M,tD> evalPoly(int deg, tC pos, DoF<M,tD> *P);

#endif

/** @} */
