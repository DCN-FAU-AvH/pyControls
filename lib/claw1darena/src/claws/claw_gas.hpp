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

/*! \file claw_gas.hpp
 *  \brief Declaration of  CLawGas for Euler gas dynamics
 */

#ifndef __GASPROBLEM_HH__
#define __GASPROBLEM_HH__

#include "../grid/grid.hpp"
#include "clawbase.hpp"
#include "bcHandler.hpp"
#include "num_source.hpp"
#include "../reconstructions/recBase.hpp"
#include "../dof/dofvector.hpp"
#include "../utils/gaussRules.hpp"
#include "../config.h"

class CLawGasBCHandler;

//! \brief CLaw class defining Euler gas dynamics equation
class CLawGas : public CLawBase<CLawGas,3>{
public:
  enum {M=3}; //!< enum for number of conserved quantities

  //! enum for the CLawGas problems
  enum PBType {PB_SIN=0, PB_SOD=1, PB_LAX=2, PB_SA=3, PB_WC=4, PB_SOD_WALL=5, PB_SLOWCONTACT=6, PB_SOD_RADIAL_2D=7, PB_HAYSAM=8 };
  //! enum for the indices of the conservative variables in t_u
  enum {RHO=0, MOM=1, EN=2};
  typedef CLawBase<CLawGas,3> t_CLawBase;
  //!< type for a set of conserved quantities defining a gas state
  typedef typename t_CLawBase::t_u t_u;

  CLawGas(tD gamma=1.4);

  static tD getGamma();

  t_u getFlMax(const t_u & u, tC, tC, tD &lMax) const;

  tD getLambdaMax(const t_u & u, tC , tC) const;

  t_u getF(const t_u & u, tC, tC) const;
  t_u getFprim(const tD & rho, const tD & vel, const tD & press) const;

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

  t_u primToCons(tD rho, tD vel, tD press) const;
  void consToPrim(t_u const & u, tD &rho, tD &vel, tD &press, tD &sound) const;

  static tD _gamma; //!< The ideal gas constant &gamma;

  void setBCFuncs(PBType pb, CLawGasBCHandler & bcH);

private:
  t_u u0sin(tC x) const;
  t_u u0sod(tC x) const;
  t_u u0lax(tC x) const;
  t_u u0shockacoustic(tC x) const;
  t_u u0slowContact(tC x) const;
  t_u u0haysam(tC x) const;

};

//! \brief CLaw class defining the Euler equation boundary treatment
class CLawGasBCHandler : public BCHandler<3>{
public:
  virtual void setLeftGhosts  (Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U);
  virtual void setRightGhosts (Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U);
private:
  bcFuncType<3> *bcFunLeft;
  bcFuncType<3> *bcFunRight;

friend void CLawGas::setBCFuncs(PBType pb , CLawGasBCHandler & bcH);
};

class CLawGasRadialSource2d : public numSource<3>{
public:
  typedef numSource<3>::t_u t_u; //!< type for a set of conserved quantities

  //! Constructor from grid and reconstruction object
  CLawGasRadialSource2d(Grid & grid, RecBase<3> & rec/*! the reconstruction object*/):
    numSource<3>(grid),
   _rec(rec)
  {}

  virtual t_u getSource(int idx, tC time) const;

  void setOrder(int order);
  int getOrder() const;

private:
  gaussRule qRule; //!< reference to the object of gaussRule
  RecBase<3> & _rec; //!< reference to the object RecBase
};

CLawGas::t_u radialSource2d(CLawGas::t_u & u, tC x, tC t);

#ifdef TD_IS_DOUBLE
/*! \brief Class for reconstruction in characteristic variables for gasdynamics
 *
 * This is a wrapper that stores internally two RecBase objects:
 * one is sized to hold the data in characteristic variables for each cell and
 * the reconstruction coefficients; the other one is sized accordingly to the
 * stencil of the RecBase and is employed during the actual computation of the
 * coefficients.
 *
 * Since the transformation between characteristic and conservative
 * variables is performed with BLAS/LAPACK, we depend on the LAPACKE library
 * and this implementation really makes sense only if tD=double!
 */
class GasCharRec: public RecBase<3>{
public:
  GasCharRec(RecBase<3>* cwiseRec);

  ~GasCharRec(){delete _sRec;}

  virtual void compute(const DofVec<3> * u, int cellStart, int cellEnd);

  virtual t_u eval(int k, tC);

  virtual const t_u * getCoeffRec(int k) const {};

  virtual void overrideCoeffRec(int k, const t_u* copy) {};

  virtual int needsGhosts() const
  {return _cwiseRec->needsGhosts();}

  virtual int getOrder() const
  {return _cwiseRec->getOrder();}

  virtual GasCharRec * duplicate(Grid & grid) const {
    std::cout<<"*** Called duplicate method of GasCharRec ***\n";
    return new GasCharRec(_cwiseRec);
  }

  virtual void setUAvgAndResize(const DofVec<3> * u){
    _uData=u;
    _vData.resize(u->size());
    _rData.resize(u->size());
    _rIPIV.resize(u->size());
    _cwiseRec->setUAvgAndResize(&_vData); //forza ridim dello storage interno alla dim di _vData, che è lungo come _uData
  }

  //! Prints info on the reconstruction
  virtual void print() const {
    std::cout << "Local Characteristic Projection, based on the ";
    _sRec->print();
  };

private:
  //this is just data storage
  std::vector<std::array<tD ,9>> _rData; //!< Matrices of right eigenvectors
  std::vector<std::array<int,3>> _rIPIV; //!< Pivot's ordering in P*rData=LU
  DofVec<3>   const * _uData; //!< pointer to the cell averages of the data in conservative variables
  DofVec<3>   _vData;         //!< cell averages of the data in characteristic variables
  RecBase<3>* _cwiseRec;      //!< RecBase where we store the rec, used to evaluate the rec

  //these are the real workers...
  Grid _sGrid;       //!< small grid with DX=1 for the data in a stencil
  DofVec<3> _sData;  //!< storage for data in the stencil of a cell
  RecBase<3>* _sRec; //!< RecBase used to compute the rec in each stencil
};

#endif

#endif

/** @} */
