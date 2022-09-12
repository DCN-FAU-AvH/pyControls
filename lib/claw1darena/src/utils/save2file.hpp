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

/** @addtogroup utils
 *
 * @{
 */

/*! \file save2file.hpp
 *  \brief Declaration of save2file function
 */

#ifndef SAVE2FILE_HH
#define SAVE2FILE_HH

#include "../dof/dofvector.hpp"
#include "../grid/grid.hpp"

template <int M>
void save2File(Grid const & G, DofVec<M> const & U, char const * filename);

#endif

/** @} */
