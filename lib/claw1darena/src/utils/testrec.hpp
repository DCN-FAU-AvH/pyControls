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

/*! \file testrec.hpp
 *  \brief Declaration of  CLawTestRec
 */

#ifndef __TESTREC_HH__
#define __TESTREC_HH__

#include "../dof/dof.hpp"
#include "gaussRules.hpp"
#include "../config.h"

//! \brief Class for the test of reconstructions
class TestRec{
public:
  enum {M=1}; //!< enum for number of conserved quantities

  enum PBType {PB_REG=0, PB_CP1=1, PB_CP2=2, PB_CP3=3, PB_CP4=4}; //!< enum for the CLawTestRec problems
  typedef DoF<M,tD> t_u;           //!< type for a set of conserved quantities

  TestRec();

  tC getCPLocation (PBType pb) const;

  void U0   (PBType pb, tC xC, t_u & u) const;
  void setU0(PBType pb, tC xC, tC dx, const gaussRule &qRule, t_u & u) const;

  PBType getPBType(const std::string& PBString) const;
  std::string getPBName(PBType pb) const;

private:

  t_u u0reg(tC x) const;
  t_u u0cp1(tC x) const;
  t_u u0cp2(tC x) const;
  t_u u0cp3(tC x) const;
  t_u u0cp4(tC x) const;

};

//! \brief Exception thrown when the name of a problem is not recognized
class UnknownProblemString{
public:
  UnknownProblemString(std::string name): _name(name) {}
  std::string _name; //!< offending string
};

//! \brief Exception thrown when a PBType is not recognized
class UnknownProblemId{
public:
  UnknownProblemId(int id): _id(id) {}
  int _id; //!< offending PBType
};

#endif

/** @} */
