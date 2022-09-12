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

#ifndef __SWEPROBLEM_TPP__
#define __SWEPROBLEM_TPP__
/*
  CLaw class defining fluxes for the shallow water equations
*/

#include <cmath>
#include <assert.h>

#include "claw_swe.hpp"

tD CLawSWE::_grav = 9.81;

//! \brief Constructor specifying the gravitational constant (g)
CLawSWE::CLawSWE(tD g)
{
  _grav = g;
}

//! \brief Returns the gravitational constant g
tD CLawSWE::getGrav()
{
  return _grav;
}

//! \brief Returns the flux and the spectral radius on the state u
typename CLawSWE::t_u CLawSWE::getFlMax(const t_u & u/*! state*/, tC, tC, tD &lMax/*! spectral radius (output)*/) const
{
  assert( u[H] >=0 );
  const tD v = u[Q]/u[H]; //TODO: desingolarizzazione della velocità
  const tD c= std::sqrt( _grav * u[H]);
  lMax = std::abs(v)+c;
  t_u f;
  f[H] = u[Q];
  f[Q] = v * u[Q] + 0.5*_grav*u[H]*u[H];
  f[Z] = 0.;
  return f;
}

//! \brief Returns the spectral radius on the state u
tD CLawSWE::getLambdaMax(const t_u & u/*! state*/, tC , tC) const
{
  assert( u[H] >=0 );
  const tD v = u[Q]/u[H]; //TODO: desingolarizzazione della velocità
  const tD c= std::sqrt( _grav * u[H]);
  return std::abs(v)+c;
}

//! \brief Returns the flux on the state u
typename CLawSWE::t_u CLawSWE::getF(const t_u & u/*! state*/, tC, tC) const
{
  assert( u[H] >=0 );
  const tD v = SWEdesingVel(u);
  t_u f;
  f[H] = u[Q];
  f[Q] = v * u[Q] + 0.5*_grav*u[H]*u[H];
  f[Z] = 0.;
  return f;
}

//! \brief Returns left domain boundary for the given PBType
tC CLawSWE::getLeftGridLimit (PBType pb) const{
  switch (pb){
  case PB_SMOOTH:
  case PB_SMOOTHFLAT:
  case PB_LAKEATREST:
  case PB_PERTLAKEATREST:
    return 0.;
  default:
    abort();
  }
}

//! \brief Returns right domain boundary for the given PBType
tC CLawSWE::getRightGridLimit(PBType pb) const
{
  switch (pb){
  case PB_SMOOTHFLAT:
  case PB_SMOOTH:
  case PB_LAKEATREST:
    return 1.;
    break;
  case PB_PERTLAKEATREST:
    return 2.;
    break;
  default:
    abort();
  }
}

//! \brief Returns final time for the given PBType
tC CLawSWE::getFinalTime(PBType pb) const
{
  switch (pb){
  case PB_SMOOTHFLAT:
  case PB_SMOOTH:
    return 0.1;
  case PB_LAKEATREST:
    return 1.0;
  case PB_PERTLAKEATREST:
    return 0.2;
  default:
    throw(UnknownProblemId(pb));  }
}

//! \brief Computes the initial data for a cell the given PBType
void CLawSWE::setU0(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const
{
  u = t_u(0.);
  t_u tmp;
  for (int k=0; k<qRule.getNNodes(); ++k){
    tC x = qRule.getNodePosition(xC,dx,k);
    switch (pb){
    case PB_SMOOTHFLAT:
      tmp =  u0smoothflat(x);
      _grav=1.0;
      break;
    case PB_SMOOTH:
      tmp =  u0smooth(x);
      _grav=1.0;
      break;
    case PB_LAKEATREST:
      tmp =  u0lakeAtRest(x);
      break;
    case PB_PERTLAKEATREST:
      tmp =  u0pertLakeAtRest(x);
      break;
    default:
      throw(UnknownProblemId(pb));
    }
    u += qRule.getWeight(k) * tmp;
  }
}

//! \brief Computes the initial data for a cell the given PBType
void CLawSWE::setUFinal(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const
{
  if (pb==PB_LAKEATREST)
    setU0(pb, xC, dx, qRule, u);
  else{
    std::cout << "Sorry, no exact solution available for problem " << getPBName(pb)
      << ". Aborting." << std::endl;
    abort();
  }
}

//! \brief Returns the the boundary condition type for the left boundary
Grid::BCType CLawSWE::getBCLeft (PBType pb) const
{
  switch (pb){
  case PB_SMOOTHFLAT:
  case PB_SMOOTH:
    return Grid::BC_PERIODIC;
  case PB_LAKEATREST:
  case PB_PERTLAKEATREST:
    return Grid::BC_FREEFLOW;
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Returns the the boundary condition type for the right boundary
Grid::BCType CLawSWE::getBCRight(PBType pb) const
{
  switch (pb){
  case PB_SMOOTHFLAT:
  case PB_SMOOTH:
    return Grid::BC_PERIODIC;
  case PB_LAKEATREST:
  case PB_PERTLAKEATREST:
    return Grid::BC_FREEFLOW;
  default:
    return getBCLeft(pb);
  }
}

//! \brief Does it have a known exact solution?
bool CLawSWE::haveExactSol(PBType pb) const{
  if (pb==PB_LAKEATREST)
    return true;
  else
    return false;
}

/*! \brief Converts string into a PBType
 *
 * Throws an UnknownProblemString(PBString) if it cannot match the
 * string with any of the official names of the problems.
 */
CLawSWE::PBType  CLawSWE::getPBType(const std::string& PBString) const
{
  if (PBString.compare("default")==0){
    return PB_SMOOTH;
  }
  else if (PBString.compare("smoothflat")==0){
    return PB_SMOOTHFLAT;
  }
  else if (PBString.compare("smooth")==0){
    return PB_SMOOTH;
  }
  else if (PBString.compare("lakeatrest")==0){
    return PB_LAKEATREST;
  }
  else if (PBString.compare("pertlakeatrest")==0){
    return PB_PERTLAKEATREST;
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
std::string CLawSWE::getPBName(PBType pb) const
{
  switch (pb){
  case PB_SMOOTHFLAT:
    return "smooth";
  case PB_SMOOTH:
    return "smooth";
  case PB_LAKEATREST:
    return "lakeatrest";
  case PB_PERTLAKEATREST:
    return "pertlakeatrest";
  default:
    throw(UnknownProblemId(pb));
  }
}

typename CLawSWE::t_u CLawSWE::u0smoothflat(tC x) const{
  t_u u;
  u[H] = 5.+std::exp(std::cos(2*M_PI*x));
  u[Q] = std::sin(std::cos(2*M_PI*x));
  u[Z] = 0.;
  return u;
}

typename CLawSWE::t_u CLawSWE::u0smooth(tC x) const{
  t_u u;
  u[H] = 5.+std::exp(std::cos(2*M_PI*x));
  u[Q] = std::sin(std::cos(2*M_PI*x));
  u[Z] = std::pow(std::sin(M_PI*x), 2);
  return u;
}

typename CLawSWE::t_u CLawSWE::u0lakeAtRest(tC x) const{
  t_u u;
  u[Z] = std::pow(std::sin(M_PI*x), 2);
  u[H] = 2. - u[Z];
  u[Q] = 0. ;
  return u;
}

typename CLawSWE::t_u CLawSWE::u0pertLakeAtRest(tC x) const{
  t_u u(0.0);
  if ( (x>1.4)&&(x<1.6) )
    u[Z] = 0.25*(1.0+std::cos(10.*M_PI*(x-0.5)));
  if ( (x>1.1)&&(x<1.2) )
    u[H] = 1.001 - u[Z];
  else
    u[H] = 1.0 - u[Z];
  u[Q] = 0. ;
  return u;
}

/*! \brief Compute v with desingularization
 *
 *  Formunlas from
 *  Kurganov, Petrova (2007)
 *  A Second-order well-balanced positivity preserving central-upwind scheme for the Saint-Venant system
 *  COMMUN. MATH. SCI., Vol. 5, No. 1, pp. 133--160
 */
tD SWEdesingVel(CLawSWE::t_u u){
  const tD tol = 1.0e-10;
  //return u[CLawSWE::H] * u[CLawSWE::Q] / (u[CLawSWE::H]*u[CLawSWE::H]+tol);

  tD h4 = u[CLawSWE::H] * u[CLawSWE::H];
  h4 = h4 * h4;
  return std::sqrt(2.0) * u[CLawSWE::H] * u[CLawSWE::Q] / std::sqrt(h4+std::max(h4,tol));
}

void CLawSWEBCHandler::setLeftGhosts(Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U){
  if (bc == Grid::BC_PERIODIC)
    setLeftPeriodicGhosts<3>(grid , U);
  else if (bc == Grid::BC_FREEFLOW)
    setLeftFreeFlowGhosts<3>(grid , U);
  else if (bc == Grid::BC_WALL){
    int k=grid.beginPhysical(); //first physical cell
    for (int i=grid.rbeginLeftGhosts(); i>grid.rendLeftGhosts() ; i--){
      U[i] = U[k];
      U[i][1] = -U[k][1]; //reflect momentum
      k++;
    }
  }
  else if (bc == Grid::BC_DIRICHLET){
    bcFuncType<3>::t_u v = (*bcFunLeft)(t);
    setLeftDirichletGhosts<3>(grid, U, v);
  } else{
    std::cout << "Error: cannot handle this boundary condition"
      << "(" << bc << ")" << std::endl;
    exit(0);
  }
}
//TODO: BC_DIRI_Q e BC_DIRI_H che impongono solo H o Q e fanno free-flow sulle altre

void CLawSWEBCHandler::setRightGhosts(Grid &grid , Grid::BCType bc, tC t, DofVec<3> &U){
  if (bc == Grid::BC_PERIODIC)
    setRightPeriodicGhosts<3>(grid , U);
  else if (bc == Grid::BC_FREEFLOW)
    setRightFreeFlowGhosts<3>(grid , U);
  else if (bc == Grid::BC_WALL){
    int k=grid.endPhysical()-1; //first physical cell
    for (int i=grid.beginRightGhosts(); i<grid.endRightGhosts() ; i++){
      U[i] = U[k];
      U[i][1] = -U[k][1]; //reflect momentum
      k--;
    }
  }
  else if (bc == Grid::BC_DIRICHLET){
    bcFuncType<3>::t_u v = (*bcFunRight)(t);
    setRightDirichletGhosts<3>(grid , U, v);
  } else{
    std::cout << "Error: cannot handle this boundary condition"
      << "(" << bc << ")" << std::endl;
    exit(0);
  }
}
//TODO: BC_DIRI_Q e BC_DIRI_H che impongono solo H o Q e fanno free-flow sulle altre

//! \brief Returns the function that represents the boundary data
void CLawSWE::setBCFuncs(CLawSWE::PBType pb, CLawSWEBCHandler & bcH){
  bcH.bcFunLeft   = NULL;
  bcH.bcFunRight  = NULL;
}

//! \brief Constructor from a base reconstruction procedure
SWEhydroRec::SWEhydroRec(RecBase<3>* cwiseRec/*! the reconstruction to be used on the characteristic variables*/):
  RecBase<3>(cwiseRec->getGrid()),  //base class for GasCharRec with full grid
  _cwiseRec(cwiseRec)               //save pointer to base rec
{};

//! \brief Computes the reconstruction coefficients in (h,q,h+z) variables
void SWEhydroRec::compute(const DofVec<3> * u, int cellStart, int cellEnd)
{
  setUAvgAndResize(u);
  const int sWidth = _cwiseRec->needsGhosts();

  for(int c=cellStart-sWidth; c<cellEnd+sWidth; ++c){
    _vData[c][0] = (*_uData)[c][CLawSWE::H];
    _vData[c][1] = (*_uData)[c][CLawSWE::Q];
    _vData[c][2] = (*_uData)[c][CLawSWE::H] + (*_uData)[c][CLawSWE::Z];
  }

  _cwiseRec->compute(&_vData, cellStart, cellEnd);
}

/*! \brief Compute the point value in conservative variables
 *
 * We first compute it in characteristic variables and then transform
 * back to conservative variables.
 *
 */
SWEhydroRec::t_u SWEhydroRec::eval(int k, tC relPos)
{
  //point-value in (h,q,h+Z) variables
  t_u rec = _cwiseRec->eval(k,relPos);
  //go back to (h,q,z)
  rec[2] -= rec[0];
  return rec;
}

void wbHydroRHS::computeRHS(t_U & U, tD time, t_U &rhs, tD &vMax) const {

    /*! - First, the boundary conditions on U are set */
  _bcHandler.setLeftGhosts (_grid , _grid.getLeftBC() , time, U);
  _bcHandler.setRightGhosts(_grid , _grid.getRightBC(), time, U);

  /*! - rhs is resized to the size of U */
  rhs.resize(U.size());
  rhs.zero();

  /*! - the reconstruction of the U data is computed */
  const int cellStart = _grid.beginPhysical()-1;
  const int cellEnd   = _grid.endPhysical()+1;
  _rec.compute(&U,cellStart,cellEnd);

  //! - compute the flux contribution of each face (using the Grid methods begin/endFaces, left/rightOfFace and isPhysical) and stack them in rhs
  vMax=0;

  for (int f = _grid.beginFaces(); f<_grid.endFaces(); ++f){
    int leftId = _grid.leftOfFace(f);
    t_u uL = _rec.evalRight(leftId);
    int rightId = _grid.rightOfFace(f);
    t_u uR = _rec.evalLeft(rightId);
    tD vmaxL, vmaxR;

    //hydrostatic well-balanced flux
    tD zFace = std::max(uL[CLawSWE::Z],uR[CLawSWE::Z]);
    tD vL = SWEdesingVel(uL);
    tD vR = SWEdesingVel(uR);
    t_u uHatL;
    tD hHatL = uL[CLawSWE::H]+uL[CLawSWE::Z]-zFace;
    uHatL[CLawSWE::H] = hHatL; uHatL[CLawSWE::Q]=hHatL*vL; uHatL[CLawSWE::Z]=zFace;
    t_u uHatR;
    tD hHatR = uR[CLawSWE::H]+uR[CLawSWE::Z]-zFace;
    uHatR[CLawSWE::H] = hHatR; uHatR[CLawSWE::Q]=hHatR*vR; uHatR[CLawSWE::Z]=zFace;
    t_u F = _numFlux.getF(uHatL, uHatR, _grid.getXFace(f), time, vmaxL, vmaxR);
    F[CLawSWE::Z] = 0.; //no flux on Z variable!

    vMax=(vmaxL>vMax? vmaxL: vMax);
    vMax=(vmaxR>vMax? vmaxR: vMax);
    if (_grid.isPhysical(leftId)){
      rhs[leftId]  -= F;
      rhs[leftId][CLawSWE::Q] += 0.5*CLawSWE::_grav*(hHatL*hHatL - uL[CLawSWE::H]*uL[CLawSWE::H]);
    }
    if (_grid.isPhysical(rightId)){
      rhs[rightId] += F;
      rhs[rightId][CLawSWE::Q]-= 0.5*CLawSWE::_grav*(hHatR*hHatR - uR[CLawSWE::H]*uR[CLawSWE::H]);
    }
  }

  //! - add the high order source discretization in the physical cells (begin/endPhysical Grid methods)
  if (_order>1)
    for (int c=_grid.beginPhysical(); c<_grid.endPhysical(); ++c){
      rhs[c] += tD(_grid.getDx()) * _source.getSource(c,time);
    }
}

int wbHydroRHS::needsGhosts() const{
  return _rec.needsGhosts() + 1;
}

int wbHydroRHS::getOrder() const{
  return _order;
}

//! Compute -g*h*Z_x using the first neighbours for the Z derivative
typename SWESourceUnb1::t_u SWESourceUnb1::getSource(int idx, tC) const{
    t_u S(0.);
    t_u recC = _rec.eval(idx  ,0.);
    t_u recL = _rec.eval(idx-1,0.);
    t_u recR = _rec.eval(idx+1,0.);

    S[CLawSWE::Q] = - CLawSWE::_grav * recC[CLawSWE::H] * (recR[CLawSWE::Z]-recL[CLawSWE::Z]) / (2.0 * _grid.getDx());
    return S;
  }

//! \brief returns the source term ..........
typename SWESourceRomberg::t_u SWESourceRomberg::getSource(int idx/*! cell index*/, tC) const{
  t_u S(0.); //inizializzo a 0

  //reconstruct at all subintervals' endpoints
  for (int q=0; q<pointLocation.size(); ++q){
    pointRec[q] = _rec.eval(idx,pointLocation[q]);
  }

  //Romberg on composite midpoint rule
  int size = _subIntervals;
  for (int l=0; l<_levels; ++l){
    tD Q = 0.; //composite midpoint of level l
    int start=0;
    while (start<_subIntervals){
      //omit -0.5*CLawSWE::_grav factor from here
      Q += (pointRec[start][CLawSWE::H]+pointRec[start+size][CLawSWE::H])
           *(pointRec[start+size][CLawSWE::Z]-pointRec[start][CLawSWE::Z]);
      start += size;
    }
    size /=2;
    //Accumulate Q in S
    S[CLawSWE::Q] += RombergCoeffs[l] * Q;
  }
  //return using -0.5*CLawSWE::_grav/dx factor
  // Note: /dx is not really needed, but we include it to be compatible with naifRHS
  return tD( - 0.5*CLawSWE::_grav / numSource<3>::_grid.getDx() ) * S;
}

/*! \brief Sets the order of the quadratur rule */
void SWESourceRomberg::setOrder(int order){
  // choose the order: smallest even larger than order
  _order= 2*((order+1)/2);

  std::cout << "Romberg of order " << _order << std::endl;

  _levels=_order/2; //number of levels
  //std::cout << " levels " << _levels << std::endl;

  // compute Romberg coeffs
  RombergCoeffs.resize(_levels);
  std::vector<tD> app;
  app.resize(RombergCoeffs.size());

  RombergCoeffs[0] = 1;
  tD den=1.;
  for (int l=1; l<_levels; ++l){
    for (int i=0; i<l ; ++i){
      app[i]=RombergCoeffs[i];
      RombergCoeffs[i]=-RombergCoeffs[i];
    }
    tD f( 1<< (2*l) );
    for (int i=0; i<l ; ++i){
      RombergCoeffs[i+1]+= f*app[i];
    }
    den *= (f-1);
  }

  //std::cout << "Romberg: ";
  for (int l=0; l<_levels; ++l){
    //std::cout << RombergCoeffs[l] << "/" << den << " " ;
    RombergCoeffs[l] /= den;
  }
  //std::cout << std::endl;

  //resize storage for inner rec points
  _subIntervals = 1<<(_levels-1); //2^(lev-1) intervals

  const int nPoints = _subIntervals + 1; //endpoints of intervals
  pointRec.resize(nPoints);
  pointLocation.resize(nPoints);
  const tD dx = 1.0 / (nPoints-1);
  for (int q=0; q<nPoints; ++q){
    pointLocation[q] = -0.5 + q*dx;
  }
}

/*! \brief Returns the order of the quadratur rule */
int SWESourceRomberg::getOrder() const{
  return _order;
}

#endif
