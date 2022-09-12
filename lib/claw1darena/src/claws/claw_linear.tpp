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

/*! \file claw_linear.tpp
 *  \brief Definition of  CLawLinear
 */

#ifndef __TRANSPORTPROBLEM_TPP__
#define __TRANSPORTPROBLEM_TPP__

/*
  CLaw class defining fluxes for linear transport with constant speed
*/

#include "../grid/grid.hpp"

#include<cmath>

#include "claw_linear.hpp"

//! \brief Default constructor (velocity=1)
CLawLinear::CLawLinear(): _vel(1.0) {}

//! \brief Constructor specifying the velocity
CLawLinear::CLawLinear(tD vel): _vel(vel) {}

//! \brief Returns the velocity
tD CLawLinear::getVel() {return _vel;}

//! \brief Returns the flux and the spectral radius on u, at x and t
typename CLawLinear::t_u CLawLinear::getFlMax
  (const t_u& u, tC, tC, tD &lMax) const
  {
    lMax = std::abs(_vel);
    return _vel * u;
  }

//! \brief Returns the spectral radius on u, at x and t
tD CLawLinear::getLambdaMax(const t_u &, tC , tC) const
{
  return std::abs(_vel);
}

//! \brief Returns the flux function on u, at x and t
typename CLawLinear::t_u CLawLinear::getF
  (const t_u& u, tC, tC) const
  {
    return _vel * u;
  }

//! \brief Returns left domain boundary for the given PBType
tC CLawLinear::getLeftGridLimit (PBType pb){
  switch (pb){
  case PB_JS:
    return -1.;
    break;
  default:
    return -0.5;
  }
}

//! \brief Returns right domain boundary for the given PBType
tC CLawLinear::getRightGridLimit (PBType pb){
  switch (pb){
  case PB_JS:
    return +1.;
    break;
  default:
    return +0.5;
  }
}

//! \brief Returns final time for the given PBType
tC CLawLinear::getFinalTime(PBType pb){
  switch (pb){
  case PB_JS:
    return 8.;
    break;
  case PB_DIRI:
    return 0.5;
    break;
  default:
    return 1.;
  }
}

//! \brief Computes the initial data for a cell the given PBType
void CLawLinear::setU0(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const{
  u = t_u(0.);
  for (int k=0; k<qRule.getNNodes(); ++k){
    switch (pb){
    case PB_SIN:
      u += qRule.getWeight(k) * u0sin(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_STEP:
      u += qRule.getWeight(k) * u0step(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_SINHF:
      u += qRule.getWeight(k) * u0sinhf(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_JS:
      u += qRule.getWeight(k) * u0JShu(qRule.getNodePosition(xC,dx,k));
      break;
    case PB_DIRI:
      u += qRule.getWeight(k) * u0zero(qRule.getNodePosition(xC,dx,k));
      break;
    default:
      throw(UnknownProblemId(pb));
    }
  }
}

//! \brief Computes the exact final data for a cell the given PBType
void CLawLinear::setUFinal(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const{
  setU0(pb, xC, dx, qRule, u);
}

//! \brief Returns the the boundary condition type for the left boundary
Grid::BCType CLawLinear::getBCLeft (PBType pb){
  if (pb==PB_DIRI)
    if (_vel>=0) return Grid::BC_DIRICHLET;
    else         return Grid::BC_FREEFLOW;
  else
    return Grid::BC_PERIODIC;
}

//! \brief Returns the the boundary condition type for the right boundary
Grid::BCType CLawLinear::getBCRight(PBType pb){
  if (pb==PB_DIRI)
    if (_vel>=0) return Grid::BC_FREEFLOW;
    else         return Grid::BC_DIRICHLET;
  else
    return Grid::BC_PERIODIC;
}

//! \brief Does it have a known exact solution?
bool CLawLinear::haveExactSol(PBType pb){
  //TODO: non è vero: si può calcolare la sol esatta anche per PB_DIRI!
  if (pb==PB_DIRI)
    return false;
  else
    return true;
}

/*! \brief Converts string into a PBType
 *
 * Throws an UnknownProblemString(PBString) if it cannot match the
 * string with any of the official names of the problems.
 */
CLawLinear::PBType CLawLinear::getPBType(const std::string& PBString) const
{
  if (PBString.compare("default")==0){
    return PB_SIN;
  }
  else if (PBString.compare("sin")==0){
    return PB_SIN;
  }
  else if (PBString.compare("step")==0){
    return PB_STEP;
  }
  else if (PBString.compare("sinhf")==0){
    return PB_SINHF;
  }
  else if (PBString.compare("jshu")==0){
    return PB_JS;
  }
  else if (PBString.compare("diri")==0){
    return PB_DIRI;
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
std::string CLawLinear::getPBName(PBType pb) const
{
  switch (pb){
  case PB_SIN:
    return("sin(2*PI*x), periodic boundary conditions");
  case PB_STEP:
    return "step, periodic boundary conditions";
  case PB_SINHF:
    return("sin(2*PI*x)+sin(30*PI*x)*exp(-80*x^2), periodic boundary conditions");
  case PB_JS:
    return("Jiang-Shu test, periodic boundary conditions");
  case PB_DIRI:
    return("test time-dependent Dirichlet boundary conditions");
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Sinusoidal initial data
typename CLawLinear::t_u CLawLinear::u0sin(tC x) const
  { return t_u(std::sin(2*M_PI*x)); }

//! \brief Double step initial data
typename CLawLinear::t_u CLawLinear::u0step(tC x) const
  { return (std::abs(x)<.25?t_u(1.0):t_u(0.)); }

//! \brief Sinusoidal high-frequency initial data
typename CLawLinear::t_u CLawLinear::u0sinhf(tC x) const
  { return t_u(std::sin(2*M_PI*x))+eWiseProduct(t_u(std::sin(30*M_PI*x)),t_u(exp(-80*pow(x,2)))); }

//! \brief Jiang-Shu initial data
typename CLawLinear::t_u CLawLinear::u0JShu(tC x) const
{
  if (x>=-0.8 & x<=-0.6)
    return ( t_u(exp(-(std::log(2) / 36 / pow(0.005,2))*pow(x-(-0.7-0.005),2))) + t_u(exp(-(std::log(2) / 36 / pow(0.005,2))*pow(x-(-0.7+0.005),2))) + tD(4) * t_u(exp(-(std::log(2) / 36 / pow(0.005,2))*pow(x+0.7,2))) ) / tD(6);
  else if (x>=-0.4 & x<=-0.2)
    return t_u(1.);
  else if (x>=0 & x<=0.2)
    return t_u(1.-std::abs(10*(x-0.1)));
  else if (x>=0.4 & x<=0.6)
    return ( t_u(std::sqrt(std::max(tD(1.-pow(10,2)*pow(x-(0.5-0.005),2)),tD(0)))) + t_u(std::sqrt(std::max(tD(1.-pow(10,2)*pow(x-(0.5+0.005),2)),tD(0)))) + tD(4.) * t_u(std::sqrt(std::max(tD(1.-pow(10,2)*pow(x-0.5,2)),tD(0.)))) ) / tD(6.);
  else
    return t_u(0.);
}

//! \brief Zero initial data
typename CLawLinear::t_u CLawLinear::u0zero(tC x) const
{ return t_u(0.); }

/*! \brief Sets the cell averages in the ghost cells on the left
 *
 * Pass in the vector of conserved quantities associated a grid and
 * on exit the left ghosts will be set to the proper values according
 * to the type of boundary condition specified.
 *
 * Only periodic and freeflow conditions are recognized.
 */
void CLawLinearBCHandler::setLeftGhosts(Grid &grid/*! grid*/, Grid::BCType bc/*! boundary condition*/, tC t, DofVec<1> &U/*! vector of conserved variables (sized according to grid)*/){
  if (bc == Grid::BC_PERIODIC)
    setLeftPeriodicGhosts<1>(grid , U);
  else if (bc == Grid::BC_FREEFLOW)
    setLeftFreeFlowGhosts<1>(grid , U);
  else if (bc == Grid::BC_DIRICHLET){
    bcFuncType<1>::t_u v = (*bcFunLeft)(t);
    setLeftDirichletGhosts<1>(grid, U, v);
  } else{
    std::cout << "Error: cannot handle this boundary condition"
      << "(" << bc << ")" << std::endl;
    exit(0);
  }
}

/*! \brief Sets the cell averages in the ghost cells on the right
 *
 * Pass in the vector of conserved quantities associated a grid and
 * on exit the right ghosts will be set to the proper values according
 * to the type of boundary condition specified.
 *
 * Only periodic and freeflow conditions are recognized.
 */
void CLawLinearBCHandler::setRightGhosts(Grid &grid/*! grid*/, Grid::BCType bc/*! boundary condition*/, tC t, DofVec<1> &U/*! vector of conserved variables (sized according to grid)*/){
  if (bc == Grid::BC_PERIODIC)
    setRightPeriodicGhosts<1>(grid , U);
  else if (bc == Grid::BC_FREEFLOW)
    setRightFreeFlowGhosts<1>(grid , U);
  else if (bc == Grid::BC_DIRICHLET){
    bcFuncType<1>::t_u v = (*bcFunRight)(t);
    setRightDirichletGhosts<1>(grid , U, v);
  } else{
    std::cout << "Error: cannot handle this boundary condition"
      << "(" << bc << ")" << std::endl;
    exit(0);
  }
}

//! sinusoidal boundary data
class bcDiriLintra: public bcFuncType<1>{
  virtual t_u operator() (tC t){ return t_u(0.5*std::sin(4.*M_PI*t)); }
};

bcDiriLintra bcDiriLintraSin;

//! \brief Returns the function that gives the boundary data for the left boundary
void CLawLinear::setBCFuncs(CLawLinear::PBType pb, CLawLinearBCHandler & bcH){
  if ((pb==CLawLinear::PB_DIRI)&&(CLawLinear::getBCLeft(pb)==Grid::BC_DIRICHLET))
    bcH.bcFunLeft  = &bcDiriLintraSin;
  else
    bcH.bcFunLeft  = NULL;
  if ((pb==CLawLinear::PB_DIRI)&&(CLawLinear::getBCRight(pb)==Grid::BC_DIRICHLET))
    bcH.bcFunRight  = &bcDiriLintraSin;
  else
      bcH.bcFunRight  = NULL;
}

#endif

/** @} */
