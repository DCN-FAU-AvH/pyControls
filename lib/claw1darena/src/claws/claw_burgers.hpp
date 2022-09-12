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

/** @addtogroup claws
 *
 * @{
 */

/*! \file claw_burgers.hpp
 *  \brief Declaration of  CLawBurgers
 */

#ifndef __BURGERSPROBLEM_HH__
#define __BURGERSPROBLEM_HH__

#include "../grid/grid.hpp"
#include "clawbase.hpp"
#include "bcHandler.hpp"
#include "../dof/dofvector.hpp"
#include "../utils/gaussRules.hpp"
#include "../config.h"

class CLawBurgersBCHandler;

//! \brief CLaw class defining the Burgers' equation
class CLawBurgers : public CLawBase<CLawBurgers,1>{
public:
  enum {M=1}; //!< enum for number of conserved quantities

  //! enum for the CLawBurgers problems
  enum PBType {PB_SIN=0,PB_STEP=1,PB_DOUBLESHOCK=2, PB_DIRI=3};
  //! enum for the indices of the conservative variables in t_u
  enum {RHO=0};

  typedef CLawBase<CLawBurgers,1> t_CLawBase;
  //!< type for a set of conserved quantities defining the state
  typedef typename t_CLawBase::t_u t_u;

  t_u getFlMax(const t_u & u, tC, tC, tD&lMax) const;

  tD getLambdaMax(const t_u & u, tC , tC) const;
  tD getVel(const t_u & u, tC , tC) const;

  t_u getF(const t_u & u, tC, tC) const;

  tC getLeftGridLimit (PBType pb) const;
  tC getRightGridLimit(PBType pb) const;
  tC getFinalTime(PBType pb) const;
  void setU0(PBType pb, tC xC, tC dx, const gaussRule &qRule, t_u & u) const;
  void setUFinal(PBType pb, tC xC, tC dx, const gaussRule &qRule, t_u & u) const;

  Grid::BCType getBCLeft (PBType pb) const;
  Grid::BCType getBCRight(PBType pb) const;

  bool haveExactSol(PBType pb) const;

  PBType getPBType(const std::string& PBString) const;
  std::string getPBName(PBType pb) const;

  void setBCFuncs(PBType pb, CLawBurgersBCHandler & bcH);

private:

  t_u u0sin(tC x) const;
  t_u u0step(tC x) const;
  t_u u0doubleshock(tC x) const;
  t_u u0zero(tC x) const;

};

//! \brief CLaw class defining the Burger's equation boundary treatment
class CLawBurgersBCHandler : public BCHandler<1>{
public:
  virtual void setLeftGhosts  (Grid &grid, Grid::BCType bc, tC t, DofVec<1> &U);
  virtual void setRightGhosts (Grid &grid, Grid::BCType bc, tC t, DofVec<1> &U);
private:
  bcFuncType<1> *bcFunLeft;
  bcFuncType<1> *bcFunRight;

friend void CLawBurgers::setBCFuncs(PBType pb , CLawBurgersBCHandler & bcH);
};

#endif

/** @} */
