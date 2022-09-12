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

/*! \file claw_linear.hpp
 *  \brief Declaration of  CLawLinear
 */

#ifndef __TRANSPORTPROBLEM_HH__
#define __TRANSPORTPROBLEM_HH__

#include "../grid/grid.hpp"
#include "clawbase.hpp"
#include "bcHandler.hpp"
#include "../dof/dofvector.hpp"
#include "../utils/gaussRules.hpp"
#include "../config.h"

class CLawLinearBCHandler;

//! \brief CLaw class defining the linear transport equation with constant speed
class CLawLinear : public CLawBase<CLawLinear,1>{
public:
  enum {M=1}; //!< enum for number of conserved quantities

  enum PBType {PB_SIN=0, PB_STEP=1, PB_SINHF=2, PB_JS=3, PB_DIRI=4}; //!< enum for the CLawLinear problems
  typedef CLawBase<CLawLinear,1> t_CLawBase;
  typedef typename t_CLawBase::t_u t_u; //!< type for a set of conserved quantities

  CLawLinear();
  CLawLinear(tD vel);

  tD getVel();

  t_u getFlMax(const t_u& u, tC, tC, tD &lMax) const;
  tD getLambdaMax(const t_u &, tC , tC) const;
  t_u getF(const t_u& u, tC, tC) const;

  tC getLeftGridLimit (PBType pb);
  tC getRightGridLimit(PBType pb);
  tC getFinalTime(PBType pb);
  void setU0(PBType pb, tC xC, tC dx, const gaussRule &qRule, t_u & u) const;
  void setUFinal(PBType pb, tC xC, tC dx, const gaussRule &qRule, t_u & u) const;
  Grid::BCType getBCLeft (PBType pb);
  Grid::BCType getBCRight(PBType pb);
  bool haveExactSol(PBType pb);

  PBType getPBType(const std::string& PBString) const;
  std::string getPBName(PBType pb) const;

  void setBCFuncs(PBType pb, CLawLinearBCHandler & bcH);

private:
  const tD _vel; //!< velocity of the transport

  t_u u0sin (tC x) const;
  t_u u0step(tC x) const;
  t_u u0sinhf(tC x) const;
  t_u u0JShu(tC x) const;
  t_u u0zero(tC x) const;

};

//! \brief CLaw class defining the linear transport equation boundary treatment
class CLawLinearBCHandler : public BCHandler<1>{
public:
  virtual void setLeftGhosts  (Grid &grid, Grid::BCType bc, tC t, DofVec<1> &U);
  virtual void setRightGhosts (Grid &grid, Grid::BCType bc, tC t, DofVec<1> &U);
private:
  bcFuncType<1> *bcFunLeft;
  bcFuncType<1> *bcFunRight;

friend void CLawLinear::setBCFuncs(PBType pb , CLawLinearBCHandler & bcH);
};

#endif

/** @} */
