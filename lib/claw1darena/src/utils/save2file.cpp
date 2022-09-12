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

//NOTA: non funziona e per ora claw1d.cpp include direttamente save2file.tpp

/** @addtogroup utils
 *
 * @{
 */

/*! \file save2file.hpp
 *  \brief Instantiation of save2file function for useful values of M
 */
#include "save2file.tpp"

template
void save2File<1>(Grid const & G, DofVec<1> const & U, char const * filename);

template
void save2File<3>(Grid const & G, DofVec<3> const & U, char const * filename);
