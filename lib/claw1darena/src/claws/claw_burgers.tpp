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

#ifndef __BURGERSPROBLEM_TPP__
#define __BURGERSPROBLEM_TPP__
/*
  CLaw class defining fluxes for the Burgers' equation
*/

#include <cmath>
#include <assert.h>

#include "claw_burgers.hpp"

//! \brief Returns the flux and the spectral radius on the state u
typename CLawBurgers::t_u CLawBurgers::getFlMax(const t_u & u/*! state*/, tC, tC, tD &lMax/*! spectral radius (output)*/) const
{
  const tD vel = u[RHO];
  lMax = std::abs(vel);
  return tD(0.5) * eWiseSquare(u);
}

//! \brief Returns the spectral radius on the state u
tD CLawBurgers::getLambdaMax(const t_u & u/*! state*/, tC , tC) const
{
  const tD vel = u[RHO];
  return std::abs(vel);
}

//! \brief Returns the velocity on the state u
tD CLawBurgers::getVel(const t_u & u/*! state*/, tC , tC) const
{
  const tD vel = u[RHO];
  return vel;
}

//! \brief Returns the flux
typename CLawBurgers::t_u CLawBurgers::getF(const t_u & u/*! state*/, tC, tC) const
{
  return tD(0.5) * eWiseSquare(u);
}

//! \brief Returns left domain boundary for the given PBType
tC CLawBurgers::getLeftGridLimit (PBType pb) const
  {return -1;}

//! \brief Returns right domain boundary for the given PBType
tC CLawBurgers::getRightGridLimit(PBType pb) const
  {return 1;}

//! \brief Returns final time for the given PBType
tC CLawBurgers::getFinalTime(PBType pb) const
{
  switch (pb){
  case PB_SIN:
    return 0.3183;
  case PB_STEP:
    return 0.5;
  case PB_DOUBLESHOCK:
    return 0.6;
  case PB_DIRI:
    return 1.0;
  default:
    throw(UnknownProblemId(pb));  }
}

//! \brief Computes the initial data for a cell the given PBType
void CLawBurgers::setU0(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const
{
  u = t_u(0.);
  t_u tmp;
  for (int k=0; k<qRule.getNNodes(); ++k){
    tC x = qRule.getNodePosition(xC,dx,k);
    switch (pb){
    case PB_SIN:
      tmp = u0sin(x);
      break;
    case PB_STEP:
      tmp = u0step(x);
      break;
    case PB_DOUBLESHOCK:
      tmp = u0doubleshock(x);
      break;
    case PB_DIRI:
      tmp = u0zero(x);
      break;
    default:
      throw(UnknownProblemId(pb));
    }
    u += qRule.getWeight(k) * tmp;
  }
}

//! \brief Computes the initial data for a cell the given PBType
void CLawBurgers::setUFinal(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const
{
  std::cout << "Sorry, no exact solution available for problem " << getPBName(pb)
    << ". Aborting." << std::endl;
  abort();
}


//! \brief Returns the the boundary condition type for the left boundary
Grid::BCType CLawBurgers::getBCLeft (PBType pb) const
{
  switch (pb){
  case PB_SIN:
  case PB_STEP:
  case PB_DOUBLESHOCK:
    return Grid::BC_PERIODIC;
  case PB_DIRI:
    return Grid::BC_DIRICHLET;
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Returns the the boundary condition type for the right boundary
Grid::BCType CLawBurgers::getBCRight(PBType pb) const
{
  switch (pb){
  case PB_SIN:
  case PB_STEP:
  case PB_DOUBLESHOCK:
    return Grid::BC_PERIODIC;
  case PB_DIRI:
    return Grid::BC_FREEFLOW;
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Does it have a known exact solution?
bool CLawBurgers::haveExactSol(PBType pb) const{
  if (pb==PB_SIN)
    return false;
  else if (pb==PB_STEP)
    return false;
  else if (pb==PB_DOUBLESHOCK)
    return false;
  else if (pb==PB_DIRI)
    return false;
  else
    return false;
}

/*! \brief Converts string into a PBType
 *
 * Throws an UnknownProblemString(PBString) if it cannot match the
 * string with any of the official names of the problems.
 */
CLawBurgers::PBType CLawBurgers::getPBType(const std::string& PBString) const
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
  else if (PBString.compare("doubleshock")==0){
    return PB_DOUBLESHOCK;
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
std::string CLawBurgers::getPBName(PBType pb) const
{
  switch (pb){
  case PB_SIN:
    return "-sin(PI*x), periodic boundary conditions";
  case PB_STEP:
    return "double step, periodic boundary conditions";
  case PB_DOUBLESHOCK:
    return "0.2-sin(PI*x)+sin(2*PI*x), periodic boundary conditions";
  case PB_DIRI:
    return "test time-dependent Dirichlet boundary conditions";
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Sinusoidal initial data
typename CLawBurgers::t_u CLawBurgers::u0sin(tC x) const
{ return t_u( - std::sin(M_PI*x) );}

//! \brief Double step initial data
typename CLawBurgers::t_u CLawBurgers::u0step(tC x) const
  { return (std::abs(x)<.25?t_u(1.0):t_u(0.)); }

//! \brief Double shock interaction test problem
typename CLawBurgers::t_u CLawBurgers::u0doubleshock(tC x) const
{ return t_u(0.2 - std::sin(M_PI*x) + std::sin(2*M_PI*x));}

//! \brief Zero initial data
typename CLawBurgers::t_u CLawBurgers::u0zero(tC x) const
{ return t_u(1.); }

void CLawBurgersBCHandler::setLeftGhosts(Grid &grid, Grid::BCType bc, tC t, DofVec<1> &U){
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

void CLawBurgersBCHandler::setRightGhosts(Grid &grid , Grid::BCType bc, tC t, DofVec<1> &U){
  if (bc == Grid::BC_PERIODIC)
    setRightPeriodicGhosts<1>(grid , U);
  else if (bc == Grid::BC_FREEFLOW)
    setRightFreeFlowGhosts<1>(grid , U);
  else if (bc == Grid::BC_DIRICHLET)
    setLeftDirichletGhosts<1>(grid , U, (*bcFunRight)(t));
  else{
    std::cout << "Error: cannot handle this boundary condition"
      << "(" << bc << ")" << std::endl;
    exit(0);
  }
}

//! cosinusoidal boundary data
class bcDiriBurger: public bcFuncType<1>{
  virtual t_u operator() (tC t){return t_u(1.0+std::sin(M_PI*t));}
};

bcDiriBurger bcDiriBurgersData;

//! \brief Returns the function that gives the boundary data for the left boundary
void CLawBurgers::setBCFuncs(CLawBurgers::PBType pb, CLawBurgersBCHandler & bcH){
  if ((pb==CLawBurgers::PB_DIRI))
    bcH.bcFunLeft  = &bcDiriBurgersData;
  else
    bcH.bcFunLeft  = NULL;
  bcH.bcFunRight  = NULL;
}

#endif
