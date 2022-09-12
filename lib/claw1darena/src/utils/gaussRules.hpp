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

/*! \file gaussRules.hpp
 *  \brief Declaration of gaussRule
 */

#include "../config.h"
#include <vector>

#ifndef GAUSSRULES_HH
#define GAUSSRULES_HH

/*! \brief Class for Gauss-Legendre quadrature rules for cell averages on an interval
 *
 * This class represents a Gauss-Legendre quadrature formula to compute
 * a cell average on an interval. The rule can be chosen by specifying the
 * number of nodes or its minimum accuracy (default=midpoint rule).
 *
 * Note that since we are dealing with cell averages,
 *     accuracy = 2 * nodes
 *
 * Usage:
 *  - instantiate an object of class gaussRule
 *  - if needed, call chooseOrder(int) or chooseNodes(int)
 *  - get the number of nodes by calling getNNodes()
 *  - for each node, get the weight by calling getWeight() and its
 * location within the cell by calling getNodePosition()
 */
class gaussRule{
public:
  gaussRule();
  void chooseOrder(int);
  void chooseNodes(int);
  int  getNNodes() const;
  int  getOrder() const;
  tC   getNodePosition(tC xC , tC H , int k) const;
  tD   getWeight(int k) const;
private:
  int _nNodes;
  std::vector<tD> _weights;
  std::vector<tC> _nodes;
};

#endif
