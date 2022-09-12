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

/** @addtogroup timeintegration
 *
 * @{
 */

/*! \file butcher.tpp
 *  \brief Declaration of ButcherTableaux
 */

#include "butcher.hpp"

#include <cassert>
#include <memory>
#include <iostream>

/*! \brief Constructor from the number of stages and the list of coefficients
 *
 * - the list is assumed to be of length S*(S+1)/2
 * - the first S*(S-1)/2 elements are interpreted as the lower
 * triangular part of A
 * - the remaining S elements are interpreted as the b coefficients
 * - the array c of abscissae is computed as row-sum of A
 */
ExplicitButcherTableaux::ExplicitButcherTableaux(int S, std::initializer_list<tD> coeffs, int order):
  _S(S), _order(order)
{
  assert(coeffs.size()==S*(S+1)/2);
  _allCoeffs = new tD[S*(S+1)/2];
  std::uninitialized_copy(coeffs.begin(), coeffs.end(), _allCoeffs);
  //Set pointers for A and b
  _A = _allCoeffs;
  _b = _allCoeffs + S*(S-1)/2;
  //store abscissae in array c
  _c = new tD[S];
  for (int i=0; i<S; ++i){
    _c[i]=0.0;
    for (int j=0;j<i;++j)
      _c[i] += getA(i,j);
  }
};

//! \brief Destructor. Deallocates storage used for the coefficients.
ExplicitButcherTableaux::~ExplicitButcherTableaux(){
  delete[] _c;
  delete[] _allCoeffs;
}

//! \brief Returns the number of stages
int ExplicitButcherTableaux::getNStages() const{
  return _S;
}

//! \brief Returns the (i,j) element of matrix A
tD ExplicitButcherTableaux::getA(int i,int j) const {
  assert(i<_S);
  assert(j<_S);
  if ( (i==0) || (j>=i))
    return 0.;
  return _A[i*(i-1)/2+j];
}

//! \brief Returns the j-th b coefficient
tD ExplicitButcherTableaux::getB(int j) const {
  assert(j<_S);
  return _b[j];
}

//! \brief Returns the i-th b coefficient
tD ExplicitButcherTableaux::getC(int i) const {
  assert(i<_S);
  return(_c[i]);
}

//! \brief Returns the order of the method
int ExplicitButcherTableaux::getOrder() const {
  return(_order);
}

//! \brief Prints out the Butcher's tableaux
void ExplicitButcherTableaux::print() const{
  for (int i=0;i<_S;++i){
    std::cout << getC(i) << " |";
    for (int j=0;j<_S;++j)
      std::cout << " " << getA(i,j);
    std::cout << std::endl;
  }
  std::cout << "---------"<<std::endl;
  std::cout << "b= ";
  for (int i=0;i<_S;++i)
    std::cout << " " << getB(i);
  std::cout << "\n---------\n";
}

/** @} */
