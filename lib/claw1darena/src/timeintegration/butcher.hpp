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

/*! \file butcher.hpp
 *  \brief Declaration of ButcherTableaux
 */

#ifndef BUTCHER_HH
#define BUTCHER_HH

#include <initializer_list>
#include "../config.h"

/*! \brief Class to store a Butcher Tableaux
 *
 * The Tableaux is assumed to be the Butcher tablueax of an
 * EXPLICIT Runge-Kutta method.
 *
 * The stages (and thus the entries in the matrix A and in the arrays
 * b and c) are numbered from 0 to S-1, where S in the number of stages.
 */
class ExplicitButcherTableaux{
public:
  ExplicitButcherTableaux(int S, std::initializer_list<tD> coeffs, int order);

  ~ExplicitButcherTableaux();

  int getNStages() const;

  tD getA(int i,int j) const;

  tD getB(int j) const;

  tD getC(int i) const;

  int getOrder() const;

  void print() const;

private:
  int _S;         //!< The number of stages
  tD *_allCoeffs; //!< Storage for the coefficients

  tD *_A;         //!< The A_ij part of the Butcher tableaux
  tD *_b;         //!< The b coefficients of the Butcher tableaux
  tD *_c;         //!< The c coefficients of the Butcher tableaux
  int _order;     //!< The order of accuracy
};

/*! \brief enum for the ERK type */
enum ERKtype {ERK_euler, ERK_heun, ERK_ssp3, ERK_runge4,ERK_rk5,ERK_rk7,ERK_rk8,ERK_rk10};

template <class S>
ExplicitButcherTableaux & getExplicitButcherTableaux(S);
//ExplicitButcherTableaux & getExplicitButcherTableaux(ERKtype type);

#endif

/** @} */
