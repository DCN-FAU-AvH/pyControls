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
 *  \file gaussRules.cpp
 *  \brief Implementation of gaussRule
 *
 * @}
 */

#include <cmath>
#include <iostream>
#include "gaussRules.hpp"

//! \brief Constructor: sets up the mid-point rule by deafult
gaussRule::gaussRule(){
  chooseNodes(1);
}

//! \brief Sets up the rule with the given number of nodes on [-0.5,0.5]
void gaussRule::chooseNodes(int nodes){
  _nNodes = nodes;
  _weights.resize(nodes);
  _nodes.resize(nodes);
  switch (nodes){
  case 1:
    _nodes  [0]=0.;
    _weights[0]=1.;
    break;
  case 2:
    _nodes  [0]= -0.5 * std::sqrt(3.)/3.0;
    _nodes  [1]= -_nodes[0];
    _weights[0]= 0.5;
    _weights[1]= 0.5;
    break;
  case 3:
    _nodes  [0]= -0.5 * std::sqrt(3./5.);
    _nodes  [1]=  0.0;
    _nodes  [2]= -_nodes[0];
    _weights[0]= 0.5 * 5./9.;
    _weights[1]= 0.5 * 8./9.;
    _weights[2]= _weights[0];
    break;
  case 4:
    _nodes  [0]= -0.5 * std::sqrt(3./7. +2./7.*std::sqrt(6./5.));
    _nodes  [1]= -0.5 * std::sqrt(3./7. -2./7.*std::sqrt(6./5.));
    _nodes  [2]= -_nodes[1];
    _nodes  [3]= -_nodes[0];
    _weights[0]= 0.5 * (18.-std::sqrt(30.))/36.;
    _weights[1]= 0.5 * (18.+std::sqrt(30.))/36.;
    _weights[2]= _weights[1];
    _weights[3]= _weights[0];
    break;
  case 5:
    _nodes  [0]= -0.5 * std::sqrt(5.+2.*std::sqrt(10./7.))/3.;
    _nodes  [1]= -0.5 * std::sqrt(5.-2.*std::sqrt(10./7.))/3.;
    _nodes  [2]= 0.0;
    _nodes  [3]= - _nodes[1];
    _nodes  [4]= - _nodes[0];
    _weights[0]= 0.5 * (322.-13.*std::sqrt(70.))/900.;
    _weights[1]= 0.5 * (322.+13.*std::sqrt(70.))/900.;
    _weights[2]= 0.5 * 128./225.;
    _weights[3]= _weights[1];
    _weights[4]= _weights[0];
    break;
  default:
    std::cout << "Gauss formula with " << nodes
      << " nodes not implemented. Aborting" << std::endl;
    abort();
  }
}

/*! \brief Sets up the rule with the minimum required accuracy on [-0.5,0.5]
 *
 * For the computation of cell averages, the Gaussian formula
 * with N nodes has accuracy 2*N
 */
void gaussRule::chooseOrder(int r){
  chooseNodes((r+1)/2);
}

//! \brief Returns the number of nodes in the current formula
int gaussRule::getNNodes() const{
  return _nNodes;
}

//! \brief Returns the order of the quadrature rule
int gaussRule::getOrder() const{
  return 2*_nNodes-1;
}

//! \brief Maps the nodes onto the given interval
tC gaussRule::getNodePosition(tC xC/*! cell center*/ , tC H /*! cell size */, int k/*! node number*/) const{
 return (xC + H * _nodes[k]);
}

//! \brief Returns a reference to the weigths for the quadrature
tD gaussRule::getWeight(int k) const{
  return  _weights[k];
}
