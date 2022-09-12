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

/** @defgroup timeintegration Time Integration
 *
 * @{
 */

/*! \file timeIntegration.hpp
 *  \brief Declaration of base class timeIntegration
 */

#ifndef TIMEINTEGRATION_HH
#define TIMEINTEGRATION_HH

#include "../config.h"
#include "../grid/grid.hpp"
#include "../dof/dofvector.hpp"
#include "../reconstructions/recBase.hpp"
#include "../numfluxes/fluxBase.hpp"
#include "../claws/num_source.hpp"
#include "../claws/bcHandler.hpp"

/*! \brief Base class for time integrators
 *
 * This class (templated on the number M of conserved variables) is a
 * base class for the time integrators of Claw1dArena.
 * It stores references to
 * - an object of class Grid (needed to get dx in timestep selection)
 *
 * and
 * - has a pure virtual method called advance that advances the
 * solution and returns the length of the timestep employed.
 *
 * In order to implement a timestepping scheme, the user should derive
 * from this class and provide an implementation of advance. The
 * following infrastructure for choosing the length of each timestep is
 * provided in this base class.
 *
 * Both fixed timestep and variable timestep are in principle allowed:
 * any implementation of advance should check the boolean variable
 * _useCFL and
 * - if _useCFL=true, advance should determine the timestep using the
 * cfl number stored in _cfl
 * - otherwise, advance should employ a timestep length of _dt
 *
 * The user will choose between the two timestepping modes by calling
 * either setCFL or setDt prior to calling advance.
 */
template<int M>
class timeIntegration{
public:
  typedef DoF<M, tD> t_u; //!< type for a set of conserved quantities
  typedef DofVec<M> t_U;  //!< type for a vector of conserved quantities on the grid

  //! \brief Constructor
  timeIntegration(Grid & grid/*! the grid*/):
    _grid(grid)
    {}

  virtual ~timeIntegration(){}

  /*! \brief Advances the solution in u0 and store the final result in u1.
   *
   * Both u0 and u1 should contain data associated to the grid stored
   * internally and set via the constructor, so that the grid methods
   * begin/endFull, begin/endPhysical etc may be used to determine which
   * are the data in u0 to be ``advanced'' in time.
   */
  virtual tC advance(t_U &u0, tC t0, t_U &u1) =0;

  //! \brief Chooses ``fixed timestep'' length advancement mode and sets dt
  void setDt(tC dt) {_dt=dt; _useCFL=false;}

  //! \brief Returns the last dt used
  tC getDt() {return _dt;}

  //! \brief Chooses ``cfl timestep'' advancement mode and sets dt
  void setCFL(tD cfl){_cfl=cfl; _useCFL=true;}

  //! \brief Returns the cfl number
  tD getCFL() {return _cfl;}

  //! \brief Sets a maximal timestep length
  void setDtMax(tC dtMax) {_dtMax=dtMax;}

  //! \brief Returns the maximal timestep length
  tC getDtMax() {return _dtMax;}

  /*! \brief Returns the number of ghosts needed by the timestepping
   * procedure. (Typically 1 o r 0)
   */
  virtual int needsGhosts() const =0;

  //! \brief Returns the theoretical order of accuracy of this timestepping procedure
  virtual int getOrder() const =0;

protected:
  tC _dt;    //!< timestep length
  tC _dtMax; //!< maximal timestep length
  tD _cfl;   //!< cfl number
  bool _useCFL;  //!< choice between ``fixed timestep'' and ``cfl timestep'' mode

  Grid& _grid;               //!< reference to the object containing the grid
};

/*! \brief Base class for right hand sides of semidiscrete methods
 *
 * This class (templated on the number M of conserved variables) is an
 * helper class for the semidiscrete time integrators of Claw1dArena.
 *
 * It has a pure virtual method called computeRHS that computes the
 * right hand side of the semidiscretized conservation/balance law and
 * returns information useful to select automatically the timestep.
 *
 * In order to implement a semidiscrete scheme, derive from this class,
 * implement the computeRHS procedure and pass an object of the derived
 * class to a timeStepping procedure like Euler<m> or explicitRungeKutta<M>
 *
 */
template<int M>
class semidiscreteRHS{
public:
  typedef DofVec<M> t_U;  //!< type for a vector of conserved quantities on the grid

  //! \brief Constructor
  semidiscreteRHS() {};

  virtual ~semidiscreteRHS() {};

  //! \brief Computes r.h.s. for the set U of conservative quantities
  virtual void computeRHS(t_U & U/*! conserved quantities*/, tD time/*! current time*/, t_U &rhs/*! (output) right hand side*/, tD &lMax/*! (output) largest eigenvalue*/ ) const =0;

  //! \brief Returns the number of ghosts needed by the timestepping procedure.
  virtual int needsGhosts()const =0;

  //! \brief Returns the theoretical order of accuracy of this timestepping procedure
  virtual int getOrder() const =0;
};

/*! \brief Naif implementation of spatial semidiscretization
 *
 * It stores references to
 * - an object of class Grid
 * - an object of some class derived from FluxBase (the numerical flux)
 * - an object of some class derived from RecBase (the  reconstrustion)
 * - an object of some class derived from BCHandler (the boundary
 * conditions)
 * - an object of some class derived from numSource (the source term)
 *
 * and
 * - computes the right hand side of the semidiscretized conservation
 * or balance law by adding up the contributions of the conservative
 * method defined by the RecBase and FluxBase objects and
 * the contribution returned by the numSource object.
 *
 * No attempt at well-balancing is done. To obtain a well-balanced scheme
 * one could select carefully the FluxBase and numSource objects passed in, or
 * write an alternative class that handles both terms toghether.
 *
 */
template<int M>
class naifRHS: public semidiscreteRHS<M>{
public:
  typedef DofVec<M> t_U;  //!< type for a vector of conserved quantities on the grid
  typedef DoF<M,tD> t_u;  //!< type for a set of conserved quantities on a cell

  //! \brief Constructor
  naifRHS(Grid & grid/*! the grid*/, RecBase<M> & rec/*! the reconstruction object*/, FluxBase<M> & numFlux/*! the numerical flux object*/, BCHandler<M>& bcHandler/*! the boundary conditions object*/, numSource<M>& source/*! the numerical source object*/):
    _grid(grid),
    _rec(rec),
    _numFlux(numFlux),
    _bcHandler(bcHandler),
    _source(source)
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
  numSource<M> & _source; //!< reference to the object providing the numerical source
  BCHandler<M> & _bcHandler; //!< reference to the object handling boundary conditions
  FluxBase<M> & _numFlux;    //!< reference to the object providing the numerical fluxes
  RecBase<M> & _rec;         //!< reference to the object providing the reconstruction
  Grid& _grid;               //!< reference to the object containing the grid

};

#endif

/** @} */
