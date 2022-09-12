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

/** @addtogroup dof
 *
 * @{
 */

/*! \file dof.cpp
 *  \brief Instantiation of DoF<M> for some value of M
 *
 *  Request instantiation of DoF<1> and DoF<3>
 */

#include "../config.h"
#include "dof.tpp"

template
  class DoF<1,tD>;
template
  DoF<1,tD> operator-(const DoF<1,tD>&);
template
  DoF<1,tD> operator+(const DoF<1,tD>&, const tD&);
template
  DoF<1,tD> operator*(const tD&, const DoF<1,tD>&);
template
  DoF<1,tD> operator/(const DoF<1,tD>&, const tD&);
template
  DoF<1,tD> operator/(const tD& , const DoF<1,tD>&);
template
  DoF<1,tD> operator-(const DoF<1,tD>&, const DoF<1,tD>&);
template
  DoF<1,tD> operator+(const DoF<1,tD>&, const DoF<1,tD>&);
template
  DoF<1,tD> eWiseMax(const DoF<1,tD>&, const DoF<1,tD>&);
template
  DoF<1,tD> eWiseAbs(const DoF<1,tD>&);
template
  DoF<1,tD> eWiseProduct(const DoF<1,tD>&, const DoF<1,tD>&);
template
  DoF<1,tD> eWiseDivide(const DoF<1,tD>&, const DoF<1,tD>&);
template
  DoF<1,tD> eWiseSquare(const DoF<1,tD>&);
template
  DoF<1,tD> eWisePow(const DoF<1,tD>&, const tD&);
template
   std::ostream& operator<<(std::ostream&, const DoF<1,tD>);
//template
//  std::istream& operator>>(std::istream& , DoF<1,tD>);

template
  class DoF<3,tD>;
template
  DoF<3,tD> operator-(const DoF<3,tD>&);
template
  DoF<3,tD> operator+(const DoF<3,tD>&, const tD&);
template
  DoF<3,tD> operator*(const tD&, const DoF<3,tD>&);
template
  DoF<3,tD> operator/(const DoF<3,tD>&, const tD&);
template
  DoF<3,tD> operator/(const tD& , const DoF<3,tD>&);
template
  DoF<3,tD> operator-(const DoF<3,tD>&, const DoF<3,tD>&);
template
  DoF<3,tD> operator+(const DoF<3,tD>&, const DoF<3,tD>&);
template
  DoF<3,tD> eWiseMax(const DoF<3,tD>&, const DoF<3,tD>&);
template
  DoF<3,tD> eWiseAbs(const DoF<3,tD>&);
template
  DoF<3,tD> eWiseProduct(const DoF<3,tD>&, const DoF<3,tD>&);
template
  DoF<3,tD> eWiseDivide(const DoF<3,tD>&, const DoF<3,tD>&);
template
  DoF<3,tD> eWiseSquare(const DoF<3,tD>&);
template
  DoF<3,tD> eWisePow(const DoF<3,tD>&, const tD&);
template
   std::ostream& operator<<(std::ostream&, const DoF<3,tD>);
//template
//  std::istream& operator>>(std::istream& , DoF<3,tD>);

/** @} */
