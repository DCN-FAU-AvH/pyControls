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

/*! \file testrec.cpp
 *  \brief Definition of TestRec
 */

#ifndef __TESTREC_CPP__
#define __TESTREC_CPP__

/*
  Class for the test of reconstructions
*/

#include<cmath>

#include "testrec.hpp"

//! \brief Default constructor
TestRec::TestRec() {}

//! \brief Computes the initial data for a point, given PBType
void TestRec::U0(PBType pb/*! PBType*/, tC x/*! location*/, t_u & u/*! cell average (output)*/) const{
  u = t_u(0.);
  switch (pb){
  case PB_REG:
    u = u0reg(x);
    break;
  case PB_CP1:
    u = u0cp1(x);
    break;
  case PB_CP2:
    u = u0cp2(x);
    break;
  case PB_CP3:
    u = u0cp3(x);
    break;
  case PB_CP4:
    u = u0cp4(x);
    break;
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Computes the initial data for a cell, given PBType
void TestRec::setU0(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const{
  u = t_u(0.);
  for (int k=0; k<qRule.getNNodes(); ++k){
    switch (pb){
    case PB_REG:
      u += qRule.getWeight(k) * u0reg(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_CP1:
      u += qRule.getWeight(k) * u0cp1(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_CP2:
      u += qRule.getWeight(k) * u0cp2(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_CP3:
      u += qRule.getWeight(k) * u0cp3(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_CP4:
      u += qRule.getWeight(k) * u0cp4(qRule.getNodePosition(xC,dx,k));
      break;
    default:
      throw(UnknownProblemId(pb));
    }
  }
}

/*! \brief Converts string into a PBType
 *
 * Throws an UnknownProblemString(PBString) if it cannot match the
 * string with any of the official names of the problems.
 */
TestRec::PBType TestRec::getPBType(const std::string& PBString) const
{
  if (PBString.compare("default")==0){
    return PB_REG;
  }
  else if (PBString.compare("reg")==0){
    return PB_REG;
  }
  else if (PBString.compare("cp1")==0){
    return PB_CP1;
  }
  else if (PBString.compare("cp2")==0){
    return PB_CP2;
  }
  else if (PBString.compare("cp3")==0){
    return PB_CP3;
  }
  else if (PBString.compare("cp4")==0){
    return PB_CP4;
  }
  else{
    std::cout << "Problema " << PBString << " non trovato!" << std::endl;
    throw(UnknownProblemString(PBString));
  }
}

/*! \brief Converts PBType into a string
 *
 * Throws an UnknownProblemId(pb) if it cannot match the
 * pb with any of the PBTypes in the enum.
 */
std::string TestRec::getPBName(PBType pb) const
{
  switch (pb){
  case PB_REG:
    return("regular function");
  case PB_CP1:
    return "critical point of order 1";
  case PB_CP2:
    return "critical point of order 2";
  case PB_CP3:
    return "critical point of order 3";
  case PB_CP4:
    return "critical point of order 4";
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Regular data: u' not zero
typename TestRec::t_u TestRec::u0reg(tC x) const
  { return t_u(std::exp(-x*x)); }

//! \brief Critical point of order 1: u'=0, u'' not zero
typename TestRec::t_u TestRec::u0cp1(tC x) const
  { return t_u(std::sin(M_PI*x-std::sin(M_PI*x)/M_PI)); }

//! \brief Critical point of order 2: u'=u''=0, u''' not zero
typename TestRec::t_u TestRec::u0cp2(tC x) const
  { return t_u( 1.0+std::pow(std::sin(M_PI*x),3) ); }

//! \brief Critical point of order 3: u'=u''=u'''=0, u'''' not zero
typename TestRec::t_u TestRec::u0cp3(tC x) const
  { return t_u( std::pow(std::cos(M_PI*x),4) ); }

//! \brief Critical point of order 4: u'=u''=u'''=u''''=0, u''''' not zero
typename TestRec::t_u TestRec::u0cp4(tC x) const
  { return t_u( 1.0+std::pow(std::sin(M_PI*x),5) ); }

/*! \brief Return location of critical point for PBType
 *
 * Throws an UnknownProblemId(pb) if it cannot match the
 * pb with any of the PBTypes in the enum.
 */
tC TestRec::getCPLocation(PBType pb) const
{
  switch (pb){
  case PB_REG:
    return 0.2;
  case PB_CP1:
    return 0.596683186911209;
  case PB_CP2:
    return 0.0;
  case PB_CP3:
    return 0.0;
  case PB_CP4:
    return 0.0;
  default:
    throw(UnknownProblemId(pb));
  }
}

#endif

/** @} */
