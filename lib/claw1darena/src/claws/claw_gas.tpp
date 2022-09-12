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

#ifndef __GASPROBLEM_TPP__
#define __GASPROBLEM_TPP__
/*
  CLaw class defining fluxes for Euler gas dynamics
*/

#include <cmath>
#include <assert.h>

#include "claw_gas.hpp"

tD CLawGas::_gamma = 1.4;

//! \brief Constructor specifying the ideal gas constant (&gamma;)
CLawGas::CLawGas(tD gamma)
{
  _gamma = gamma;
}

//! \brief Returns the ideal gas constant &gamma;
tD CLawGas::getGamma()
{
  return _gamma;
}

//! \brief Returns the flux and the spectral radius on the state u
typename CLawGas::t_u CLawGas::getFlMax(const t_u & u/*! state*/, tC, tC, tD &lMax/*! spectral radius (output)*/) const
{
  assert( u[RHO] >=0 );
  const tD v = u[MOM]/u[RHO];
  const tD p = (_gamma-1.0) * (u[EN] - 0.5*u[MOM]*v);
  assert ( p >=0 );
  const tD c= std::sqrt( _gamma * p / u[RHO]);
  lMax = std::abs(v)+c;
  t_u f;
  f[RHO] = u[MOM];
  f[MOM] = v * u[MOM] + p;
  f[EN]  = v * (u[EN] + p);
  return f;
}

//! \brief Returns the spectral radius on the state u
tD CLawGas::getLambdaMax(const t_u & u/*! state*/, tC , tC) const
{
  assert( u[RHO] >=0 );
  const tD v = u[MOM]/u[RHO];
  const tD p = (_gamma-1.0) * (u[EN] - 0.5*u[MOM]*v);
  assert ( p >=0 );
  const tD c= std::sqrt( _gamma * p / u[RHO]);
  return (std::abs(v)+c);
}

//! \brief Returns the flux on the state u
typename CLawGas::t_u CLawGas::getF(const t_u & u/*! state*/, tC, tC) const
{
  assert( u[RHO] >=0 );
  const tD v = u[MOM]/u[RHO];
  const tD p = (_gamma-1.0) * (u[EN] - 0.5*u[MOM]*v);
  assert ( p >=0 );
  t_u f;
  f[RHO] = u[MOM];
  f[MOM] = v * u[MOM] + p;
  f[EN]  = v * (u[EN] + p);
  return f;
}

//! \brief Returns the flux on a state given in primitive variables
typename CLawGas::t_u CLawGas::getFprim(const tD & rho/*! density*/, const tD & vel/*! velocity*/, const tD & press/*! pressure*/) const
{
  t_u f;
  f[RHO] = rho * vel;
  f[MOM] = rho * vel * vel + press;
  f[EN]  = vel * (0.5*rho*vel*vel + _gamma/(_gamma-1.) * press);
  return f;
}

//! \brief Returns left domain boundary for the given PBType
tC CLawGas::getLeftGridLimit (PBType pb) const{
  if (pb==PB_SLOWCONTACT)
    return -5.;
  else if (pb==PB_SA)
    return -5.;
  else
    return 0.;
}

//! \brief Returns right domain boundary for the given PBType
tC CLawGas::getRightGridLimit(PBType pb) const
{
  if (pb==PB_SLOWCONTACT)
    return 5.;
  else if (pb==PB_SA)
    return 5.;
  else
    return 1.;
}

//! \brief Returns final time for the given PBType
tC CLawGas::getFinalTime(PBType pb) const
{
  switch (pb){
  case PB_SIN:
    return 1.0;
  case PB_SOD:
    return 0.2;
  case PB_LAX:
    return 0.15;
  case PB_SA:
    return 1.8;
  case PB_WC:
    return 0.04;
  case PB_SOD_WALL:
    return 0.5;
  case PB_SLOWCONTACT:
    return 10.0;
  case PB_SOD_RADIAL_2D:
    return 0.25;
  case PB_HAYSAM:
    return 2.0;
  default:
    throw(UnknownProblemId(pb));  }
}

//! \brief Computes the initial data for a cell the given PBType
void CLawGas::setU0(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const
{
  u = t_u(0.);
  t_u tmp;
  for (int k=0; k<qRule.getNNodes(); ++k){
    tC x = qRule.getNodePosition(xC,dx,k);
    switch (pb){
    case PB_SIN:
      tmp =  u0sin(x);
      break;
    case PB_SOD:
    case PB_SOD_WALL:
    case PB_SOD_RADIAL_2D:
      tmp = u0sod(x);
      break;
    case PB_LAX:
      tmp = u0lax(x);
      break;
    case PB_SA:
      tmp = u0shockacoustic(x);
      break;
    case PB_SLOWCONTACT:
      tmp = u0slowContact(x);
      break;
    case PB_HAYSAM:
      tmp = u0haysam(x);
      break;
    default:
      throw(UnknownProblemId(pb));
    }
    u += qRule.getWeight(k) * tmp;
  }
}

//! \brief Computes the initial data for a cell the given PBType
void CLawGas::setUFinal(PBType pb/*! PBType*/, tC xC/*! cell center*/, tC dx/*! cell size*/, const gaussRule &qRule/*! quadrature rule*/, t_u & u/*! cell average (output)*/) const
{
  if (pb==PB_SIN)
    setU0(pb, xC, dx, qRule, u);
  else{
    std::cout << "Sorry, no exact solution available for problem " << getPBName(pb)
      << ". Aborting." << std::endl;
    abort();
  }
}


//! \brief Returns the the boundary condition type for the left boundary
Grid::BCType CLawGas::getBCLeft (PBType pb) const
{
  switch (pb){
  case PB_SIN:
    return Grid::BC_PERIODIC;
  case PB_SOD:
  case PB_LAX:
  case PB_SA:
  case PB_SLOWCONTACT:
  case PB_SOD_RADIAL_2D:
    return Grid::BC_FREEFLOW;
  case PB_WC:
  case PB_SOD_WALL:
    return Grid::BC_WALL;
  case PB_HAYSAM:
    return Grid::BC_DIRICHLET;
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Returns the the boundary condition type for the right boundary
Grid::BCType CLawGas::getBCRight(PBType pb) const
{
  switch (pb){
  case PB_HAYSAM:
    return Grid::BC_WALL;
  default:
    return getBCLeft(pb);
  }
}

//! \brief Does it have a known exact solution?
bool CLawGas::haveExactSol(PBType pb) const{
  if (pb==PB_SIN)
    return true;
  else
    return false;
}

/*! \brief Converts string into a PBType
 *
 * Throws an UnknownProblemString(PBString) if it cannot match the
 * string with any of the official names of the problems.
 */
CLawGas::PBType  CLawGas::getPBType(const std::string& PBString) const
{
  if (PBString.compare("default")==0){
    return PB_SIN;
  }
  else if (PBString.compare("sin")==0){
    return PB_SIN;
  }
  else if (PBString.compare("sod")==0){
    return PB_SOD;
  }
  else if (PBString.compare("lax")==0){
    return PB_LAX;
  }
  else if (PBString.compare("shockacoustic")==0){
    return PB_SA;
  }
  else if (PBString.compare("sodwall")==0){
    return PB_SOD_WALL;
  }
  else if (PBString.compare("slowcontact")==0){
    return PB_SLOWCONTACT;
  }
  else if (PBString.compare("sodradial2d")==0){
    return PB_SOD_RADIAL_2D;
  }
  else if (PBString.compare("haysam")==0){
    return PB_HAYSAM;
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
std::string CLawGas::getPBName(PBType pb) const
{
  switch (pb){
  case PB_SIN:
    return "sinusoidal";
  case PB_SOD:
    return "Sod shock tube";
  case PB_LAX:
    return "Lax shock tube";
  case PB_SA:
    return "Shock-acoustic interaction";
  case PB_SOD_WALL:
    return "Sod shock tube with wall boundary conditions";
  case PB_SLOWCONTACT:
    return "Slow, right-moving contact discontinuity";
  case PB_SOD_RADIAL_2D:
    return "Sod's explosion problem in spherical coordinates";
  case PB_HAYSAM:
    return "Haysam's problem";
  default:
    throw(UnknownProblemId(pb));
  }
}

//! \brief Converts primitive variable into a set of conservative variables stored in a t_u
typename CLawGas::t_u CLawGas::primToCons(tD rho/*! density*/, tD vel/*! velocity*/, tD press/*! pressure */) const
{
  t_u u;
  u[RHO] = rho;
  u[MOM] = rho*vel;
  u[EN]  = 0.5*rho*vel*vel + press/(_gamma-1.);
  return u;
}

//! \brief Converts a set of conservative variables stored in a t_u into primitive variables
void CLawGas::consToPrim(t_u const & u/*! conservative variables*/, tD &rho/*! density (output)*/, tD &vel/*! velocity (output)*/, tD &press/*! pressure (output)*/, tD &sound/*! sound speed (output) */) const
{
  rho = u[RHO];
  vel = u[MOM]/rho;
  press = (_gamma-1.) * (u[EN] - 0.5*rho*vel*vel);
  sound = std::sqrt(_gamma*press/rho);
}

typename CLawGas::t_u CLawGas::u0sin(tC x) const
{ return primToCons(1+0.5*std::sin(2*M_PI*x), 1.0, 1.0);}

typename CLawGas::t_u CLawGas::u0sod(tC x) const
{
  if (x<0.5)
    return primToCons(1.0  , 0.0, 1.0);
  else
    return primToCons(0.125, 0.0, 0.1);
}

typename CLawGas::t_u CLawGas::u0lax(tC x) const
{
  if (x<0.5)
    return primToCons(0.445, 0.698, 3.528);
  else
    return primToCons(0.5  , 0.0  , 0.571);
}

typename CLawGas::t_u CLawGas::u0shockacoustic(tC x) const
{
  if (x<-4.0)
    return primToCons(3.857143, 2.629369, 10.333333);
  else
    return primToCons(1+0.2*std::sin(5*x), 0.0, 1.0);
  }

typename CLawGas::t_u CLawGas::u0haysam(tC x) const
{
  return primToCons(1.0  , 0.0, 1.0);
}

typename CLawGas::t_u CLawGas::u0slowContact(tC x) const
{
  if (x<0.0)
    return primToCons(2.0, 0.1, 1.0);
   else
    return primToCons(1.0, 0.1, 1.0);
}

void CLawGasBCHandler::setLeftGhosts(Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U){
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

void CLawGasBCHandler::setRightGhosts(Grid &grid , Grid::BCType bc, tC t, DofVec<3> &U){
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

//! Haysam's problem boundary data
class bcDiriHaysam: public bcFuncType<3>{
  virtual t_u operator() (tC t){
    const tD vel=(t<0.5?tD(0.01)*std::pow(std::sin(tD(2.*M_PI)*t) , tD(3.0)):tD(0.0));
    const tD rho=1.0;
    const tD press=1.0;
    t_u u;
    u[CLawGas::RHO] = rho;
    u[CLawGas::MOM] = rho*vel;
    u[CLawGas::EN]  = 0.5*rho*vel*vel + press/(CLawGas::_gamma-1.);
    return u;
  }
};

bcDiriHaysam bcDiriHaysamData;

//! \brief Returns the function that represents the boundary data
void CLawGas::setBCFuncs(CLawGas::PBType pb, CLawGasBCHandler & bcH){
  if (pb==CLawGas::PB_HAYSAM)
    bcH.bcFunLeft  = &bcDiriHaysamData;
  else
    bcH.bcFunLeft  = NULL;
  bcH.bcFunRight  = NULL;
}

//! \brief returns the source term for the radial 2d Euler equations
typename CLawGasRadialSource2d::t_u CLawGasRadialSource2d::getSource(int idx/*! cell index*/, tC time) const{
  t_u S;
  const tD xc = numSource<3>::_grid.getXCenter(idx);
  const tD dx = numSource<3>::_grid.getDx();
  for (int k=0; k<qRule.getNNodes(); ++k){
    tC xq = qRule.getNodePosition(xc,dx,k);
    t_u urec = _rec.eval(idx,(xq-xc)/dx);
    S += qRule.getWeight(k) * radialSource2d(urec,xq,time);
  }
  return S;
}

//! \brief returns the source term for the radial 2d Euler equations
typename CLawGas::t_u radialSource2d(CLawGas::t_u & u/*! state*/, tC x/*! position*/, tC t/*! time*/)
{
  assert( u[CLawGas::RHO] >=0 );
  const tD v = u[CLawGas::MOM]/u[CLawGas::RHO];
  const tD p = (CLawGas::_gamma-1.0) * (u[CLawGas::EN] - 0.5*u[CLawGas::MOM]*v);
  assert ( p >=0 );
  CLawGas::t_u S;
  tD alpha=1;
  S[CLawGas::RHO] = - alpha/x * u[CLawGas::MOM];
  S[CLawGas::MOM] = - alpha/x * u[CLawGas::MOM] * u[CLawGas::MOM]/u[CLawGas::RHO];
  S[CLawGas::EN]  = - alpha/x * u[CLawGas::MOM]/u[CLawGas::RHO] * (u[CLawGas::EN] + p);
  return S;
}

/*! \brief Sets the order of the quadratur rule */
void CLawGasRadialSource2d::setOrder(int order){
  qRule.chooseOrder(order);
}

/*! \brief Returns the order of the quadratur rule */
int CLawGasRadialSource2d::getOrder() const{
  return qRule.getOrder();
}

#if defined(GASDIN) && defined(TD_IS_DOUBLE)
// This code relies on lapacke and cblas which are linked only into gasdynamics.
//   However this file is included by claw1d.ccp also when compiling for other
// conservation laws, which in turns causes a linking error. For now we only
// compile the code for characteristic reconstructions only when GASDIN is
// defined, which is a horrible but effective hack...
// The proper solution is to fix the compilation of the numerical fluxes,
// maybe moving their implementation in the claw_***.tpp files.
//
// TODO!
#include "lapacke.h"
#include "cblas.h"

//! \brief Constructor from a base reconstruction procedure
GasCharRec::GasCharRec(RecBase<3>* cwiseRec/*! the reconstruction to be used on the characteristic variables*/):
  RecBase<3>(cwiseRec->getGrid()),  //base class for GasCharRec with full grid
  _cwiseRec(cwiseRec),              //save pointer to base rec
  _sGrid(0., 1., 1, cwiseRec->needsGhosts()) //reduced grid for stencil data: 1 cell + ghosts
{
  _sRec = _cwiseRec->duplicate(_sGrid);
  _sData.resize(_sGrid.getFullSize());
};

/*! \brief Computes the reconstruction coefficients in characteristic variables
 *
 * Leaves in _rData and _rIPIV the needed data to transform back to
 * conservativa variables after evaluating the reconstruction.
 */
void GasCharRec::compute(const DofVec<3> * u, int cellStart, int cellEnd)
{
  const tD & gamma = CLawGas::_gamma;
  setUAvgAndResize(u);
  const int sWidth = _cwiseRec->needsGhosts();
  _sRec->setUAvgAndResize(&_sData);
  for(int i=cellStart; i<cellEnd; ++i){
    //copy stencil data in temp vector _sData
    for (int k=-sWidth; k<=sWidth; ++k)
      _sData[k+sWidth] = (*u)[i+k];
    //index of central cell in the local grid
    const int myCell = _sGrid.beginPhysical();

    //! The code is based on the fact that in the DofVec _sData
    //! _sData[0].data() points to a contiguous memory area with all
    //! the data stored as
    //! _sData[0][0] _sData[0][1] _sData[0][2] _sData[1][0] ...
    //! i.e. in col-major format

    // Matrix to tranform data (col-major format)
    const t_u & uRef = (*u)[i];
    const tD v= uRef[CLawGas::MOM]/uRef[CLawGas::RHO];
    const tD p=(gamma-1.)*(uRef[CLawGas::EN]-0.5*uRef[CLawGas::RHO]*v*v);
    const tD c=std::sqrt(gamma*p/uRef[CLawGas::RHO]);
    const tD H=(uRef[CLawGas::EN]+p)/uRef[CLawGas::RHO];
    _rData[i] = {1,v-c,H-v*c,  1,v,0.5*v*v, 1,v+c,H+v*c};

    //Transform stencil data
    const int NRHS = _sGrid.getFullSize();
    int INFO = LAPACKE_dgesv_work(LAPACK_COL_MAJOR, 3, NRHS, _rData[i].data(),  3, _rIPIV[i].data(), _sData[0].data(), 3);
    //= 0:  successful exit
    //< 0:  if INFO = -i, the i-th argument had an illegal value
    //> 0:  if INFO = i, U(i,i) in the factorization is exactly zero
    if (INFO){
      std::cout<<"Error " << INFO << " in LAPACKE_degsv_work in characteristic projection for cell " << i << std::endl;
      abort();
    }
    //! As a result of calling LAPACKE_dgesv_work to compute the characteristic projection,
    //! upon exit from this function, _rData will contain the PA=LU decomposition of the matrix
    //! of right eigenvectors and _rIPIV the order of pivotal elements employed in the decomposition.

    //Store avrage of central cell in characteristic variables
    _vData[i] = _sData[myCell];
    //Compute the reconstruction coefficients in characterictic variables
    _sRec->compute(&_sData,myCell,myCell+1);
    //Store the coefficients
    _cwiseRec->overrideCoeffRec(i,_sRec->getCoeffRec(myCell)); //copia i coeff della ricostruzione
  }
}

/*! \brief Compute the point value in conservative variables
 *
 * We first compute it in characteristic variables and then transform
 * back to conservative variables.
 *
 */
GasCharRec::t_u GasCharRec::eval(int k, tC relPos)
{
  //point-value in characteristic variables
  t_u rec = _cwiseRec->eval(k,relPos);

  //transform back to conservative variables
  //! The product R*rec is computed as P^{-1}*L*U*rec, using data left
  //! in _rData and _rIPIV by the compute method
  cblas_dtrmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 3, _rData[k].data(), 3, rec.data(), 1); //U*x
  cblas_dtrmv(CblasColMajor, CblasLower, CblasNoTrans, CblasUnit   , 3, _rData[k].data(), 3, rec.data(), 1); //L*x
  LAPACKE_dlaswp_work(LAPACK_COL_MAJOR, 1, rec.data(), 3, 1, 3, _rIPIV[k].data(), -1 );                      //P^{-1}*x

  return rec;
}

#endif

#endif
