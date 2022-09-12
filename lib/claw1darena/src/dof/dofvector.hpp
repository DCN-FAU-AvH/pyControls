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

/*! \file dofvector.hpp
 *  \brief Declaration of DoFVec<M>
 */

#ifndef __DOFVECTOR_HH
#define __DOFVECTOR_HH

#include <vector>
#include "dof.hpp"

/*! \brief A class to store the cell averages on the grid
 *
 * A std::vector<DoF<M>> that exposes some methods of
 * a std::vector and with a zero() method to sero out the
 * entire vector.
*/

template <int M>
class DofVec{
private:
  /*! \brief type for DoF*/
  typedef DoF<M, tD> t_u;
  /*! \brief type for internal data storage*/
  typedef typename std::vector<t_u> t_U;
  /*! \brief internal data storage*/
  t_U U;
public:
  typedef typename t_U::size_type size_t;
  //typedef t_u value_t;
  /*! \brief Default constructor */
  DofVec():U(){}
  /*! \brief Constructs a DofVec of a given size */
  DofVec(int n):U(n){}
  /*! \brief Access to element */
  t_u & operator[](int i)
    {return U[i];}
  /*! \brief Const access to element */
  const t_u & operator[](int i) const
    {return U[i];}
  /*! \brief Resize the DofVec */
  void resize(int i){U.resize(i);}
  /*! \brief Reserve space for the DofVec */
  void reserve(int i){U.reserve(i);}
  /*! \brief Returns the size the DofVec */
  size_t size() const {return U.size();}
  /*! \brief Swap data with another DofVec */
  void swap(DofVec U2){U.swap(U2.U);}
  /*! \brief Zero out all elements of the DofVec */
  void zero() {std::fill(U.begin(),U.end(),t_u(0.));}
};

#endif

/** @} */
