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

/*! \file claw_swe.hpp
 *  \brief Declaration of  CLawSWE for the shallow water equations
 *
 * This is defined as a system of conservation laws with 3 equations,
 * where the last one is dz/dt=0 so that all vectors of conserved
 * quantities, reconstructions, etc have size 3. This wastes a bit of
 * memory since of course z never changes, but saves a few
 * headaches...
 */

#ifndef __SWEPROBLEM_HH__
#define __SWEPROBLEM_HH__

#include "../config.h"
#include "../grid/grid.hpp"
#include "../grid/gridutils.hpp"
#include "clawbase.hpp"
#include "bcHandler.hpp"
#include "num_source.hpp"
#include "../reconstructions/recBase.hpp"
#include "../dof/dofvector.hpp"
#include "../timeintegration/timeIntegration.hpp"
#include "../utils/gaussRules.hpp"

class CLawSWEBCHandler;

//! \brief CLaw class defining the shallow water equation
class CLawSWE : public CLawBase<CLawSWE,3>{
public:
  enum {M=3}; //!< enum for number of conserved quantities

  //! enum for the CLawSWE problems
  enum PBType {PB_SMOOTH=0, PB_SMOOTHFLAT=1, PB_LAKEATREST=2, PB_PERTLAKEATREST=3};
  //! enum for the indices of the conservative variables in t_u
  enum {H=0, Q=1, Z=2};
  typedef CLawBase<CLawSWE,3> t_CLawBase;
  //!< type for a set of conserved quantities defining a SWE state
  typedef typename t_CLawBase::t_u t_u;

  CLawSWE(tD grav=9.81);

  static tD getGrav();

  t_u getFlMax(const t_u & u, tC, tC, tD &lMax) const;

  tD getLambdaMax(const t_u & u, tC , tC) const;

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

  static tD _grav; //!< The gravitational constant g

  void setBCFuncs(PBType pb, CLawSWEBCHandler & bcH);

private:
  t_u u0smooth(tC x) const;
  t_u u0smoothflat(tC x) const;
  t_u u0lakeAtRest(tC x) const;
  t_u u0pertLakeAtRest(tC x) const;

};

/*! \brief Class for hydrostatic reconstruction for Saint-Venant equations
 *
 */
class SWEhydroRec: public RecBase<3>{
public:
  SWEhydroRec(RecBase<3>* cwiseRec);

  virtual void compute(const DofVec<3> * u, int cellStart, int cellEnd);

  virtual t_u eval(int k, tC);

  virtual const t_u * getCoeffRec(int k) const {};

  virtual void overrideCoeffRec(int k, const t_u* copy) {};

  virtual int needsGhosts() const
  {return _cwiseRec->needsGhosts();}

  virtual int getOrder() const
  {return _cwiseRec->getOrder();}

  virtual SWEhydroRec * duplicate(Grid & grid) const {
    std::cout<<"*** Called duplicate method of SWEhydroRec ***\n";
    return new SWEhydroRec(_cwiseRec);
  }

  virtual void setUAvgAndResize(const DofVec<3> * u){
    _uData=u;
    _vData.resize(u->size());
    _cwiseRec->setUAvgAndResize(&_vData); //forza ridim dello storage interno alla dim di _vData, che è lungo come _uData
  }

  //! Prints info on the reconstruction
  virtual void print() const {
    std::cout << "SWE hydrostatic, based on the ";
    _cwiseRec->print();
  };

private:
  //this is just data storage
  DofVec<3>   const * _uData; //!< pointer to the cell averages of the data in conservative variables
  DofVec<3>   _vData;         //!< cell averages of the data in characteristic variables
  RecBase<3>* _cwiseRec;      //!< RecBase where we store the rec, used to evaluate the rec
};

/*! \brief SWE source unbalanced central quadrature for h*Z_x
 *
 * The term Z_x is computed with central differencing using the values
 * at the cell centers of the first neightbours.
 *
 * This class is really thought for computing a 1st order accurate source
 * when using constant reconstructions. When using any nontrivial
 * reconstruction, the SWESourceRomberg class should be used.
 */
class SWESourceUnb1 : public numSource<3>{
public:
  typedef numSource<3>::t_u t_u; //!< type for a set of conserved quantities

  //! Constructor from grid and reconstruction object
  SWESourceUnb1(Grid & grid, RecBase<3> & rec/*! the reconstruction object*/):
    numSource<3>(grid),
    _rec(rec)
    {}

  virtual t_u getSource(int idx, tC time) const;

  int getOrder() const
    {return 1;}

private:
  RecBase<3> & _rec; //!< reference to the object RecBase
};

/*! \brief SWE source quadrature with Romberg extrapolation of midpoint rule for h*Z_x
 *
 * This implements ideas from
 * Noelle, Pankratz, Puppo, Natvig
 * Well-balanced finite volume schemes of arbitrary order of accuracy for shallow water flows
 * Journal of Computational Physics 213 (2006) 474--499
 * doi:10.1016/j.jcp.2005.08.019
 */
class SWESourceRomberg : public numSource<3>{
public:
  typedef numSource<3>::t_u t_u; //!< type for a set of conserved quantities

  //! Constructor from grid and reconstruction object
  SWESourceRomberg(Grid & grid, RecBase<3> & rec/*! the reconstruction object*/):
    numSource<3>(grid),
   _rec(rec)
  {setOrder(rec.getOrder());}

  virtual t_u getSource(int idx, tC time) const;

  int getOrder() const;

private:
  RecBase<3> & _rec; //!< reference to the object RecBase
  int _order;
  int _levels;
  int _subIntervals;
  std::vector<tD>  RombergCoeffs;
  mutable std::vector<t_u> pointRec;
  std::vector<tD>  pointLocation;

  void setOrder(int order);
};

/*! \brief Well-balanced hydrostatic method
 *
 * It stores references to
 * - an object of class Grid
 * - an object of some class derived from FluxBase (the numerical flux)
 * - an object of some class derived from RecBase (the  reconstruction)
 * - an object of some class derived from BCHandler (the boundary
 * conditions)
 *
 * and
 * - computes the right hand side of the semidiscretized conservation
 * or balance law by performing the reconstruction with the
 * given RecBase object, the given FluxBase on the modified
 * boundary data and adding up the contribution of a high order
 * well-balanced quadrature based on Romberg extrapolation of the
 * midpoint rule for h*Z_x.
 *
 * Note: this will give a well-balanced method on the lake-at-rest
 * equilibrium if the RecBase object computes a globally continuous
 * reconstruction on such steady state. For example a SWEhydroRec
 * object will do the job.
 *
 * This implements ideas from
 * Audusse, Bouchut, Bristeau, Klein, Perthame
 * A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows
 * SIAM J. SCI. COMPUT. Vol. 25 (2004), No. 6, pp. 2050--2065
 * doi: 10.1137/S1064827503431090
 */
class wbHydroRHS: public semidiscreteRHS<3>{
public:
  typedef DofVec<3> t_U;  //!< type for a vector of conserved quantities on the grid
  typedef DoF<3,tD> t_u;  //!< type for a set of conserved quantities on a cell

  //! \brief Constructor
  wbHydroRHS(Grid & grid/*! the grid*/, RecBase<3> & rec/*! the reconstruction object*/, FluxBase<3> & numFlux/*! the numerical flux object*/, BCHandler<3>& bcHandler/*! the boundary conditions object*/):
    _grid(grid),
    _rec(rec),
    _numFlux(numFlux),
    _bcHandler(bcHandler),
    _order(rec.getOrder()),
    _source(_grid,_rec)
    {}

  //! \brief Computes r.h.s. for the set U of conservative quantities
  virtual void computeRHS(t_U & U, tD time, t_U &rhs, tD &lMax) const;

  /*! \brief Returns the number of ghosts needed by the timestepping
   * procedure. (Typically 1 o r 0)
   */
  virtual int needsGhosts() const;

  //! \brief Returns the theoretical order of accuracy of this timestepping procedure
  virtual int getOrder() const;

protected:
  BCHandler<3> & _bcHandler; //!< reference to the object handling boundary conditions
  FluxBase<3> & _numFlux;    //!< reference to the object providing the numerical fluxes
  RecBase<3> & _rec;         //!< reference to the object providing the reconstruction
  Grid& _grid;               //!< reference to the object containing the grid

private:
  //mutable DofVec<3>  _vData;         //!< cell averages of the data in (h,q,h+Z) variables
  const int _order;
  SWESourceRomberg _source;
};

//! \brief CLaw class defining the shallow water equation boundary treatment
class CLawSWEBCHandler : public BCHandler<3>{
public:
  virtual void setLeftGhosts  (Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U);
  virtual void setRightGhosts (Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U);
private:
  bcFuncType<3> *bcFunLeft;
  bcFuncType<3> *bcFunRight;

friend void CLawSWE::setBCFuncs(PBType pb , CLawSWEBCHandler & bcH);
};

tD SWEdesingVel(CLawSWE::t_u u);

#endif

/** @} */
