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

/*! \file udivdiff.hpp
 *  \brief Declaration of UDivDiff
 */

#ifndef UDIVDIFF_HH
#define UDIVDIFF_HH

#include "../config.h"
#include "../dof/dof.hpp"
#include "../dof/dofvector.hpp"
#include <vector>
#include <array>

/*! \brief Class to compute, store and use (Un)Divided Differences
 *
 * It is templated on the number of conserved quantities M.
 *
 * The method compute takes in input a vector of conserved quantities
 * and fills in the matrix of (un)divided differences up to the
 * specified order.
 *
 * IMPORTANT: A pointer to the input vector is stored internally for
 * later use. Such data should not change between the call to compute
 * and the call to functions like InterPoly that makes use of the input
 * cell averages!
 *
 * Other methods makes clever use of the computed table to return the
 * coefficients of polynomials interpolating the cell averages passed
 * to the compute method.
 */
template <int M>
class UDivDiff{
public:
  typedef DoF<M, tD> t_u; //!< type for a set of conserved quantities
  UDivDiff();
  void compute(DofVec<M> const & u,int maxD);
  t_u getUDD(int cell, int deg);
  void InterpPoly(int cell, int deg, int shift,t_u* P);

private:
  int _maxD;
  std::vector<std::vector<t_u>> _udd; //!< vector to store the (un)divided differences
  const DofVec<M>* _u; //!< pointer to the cell averages last passed to compute()
};

#endif

/** @} */
