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

/** @defgroup reconstructions Reconstructions
 *
 * @{
 */

/*! \file clawbase.hpp
 *  \brief Declaration of CLawBase
 */

#ifndef RECBASE_HH
#define RECBASE_HH

#include "../config.h"
#include "../dof/dofvector.hpp"
#include "../grid/grid.hpp"

/*! \brief Base class for reconstructions
 *
 * Pure virtual class, templated on the number of conservative variables.
 *
 * In order to define a new reconstruction, derive a class from this one
 * and implement the compute, eval, evalL, evalR methods.
 */

template <int M>
class RecBase{
public:
  typedef DoF<M, tD> t_u; //!< type for a set of conserved quantities

  /*! \brief Default constructor
   *
   * Stores a reference to the grid in the internal variable.
   */
  RecBase(Grid & grid):
    _grid(grid)
    {}

  virtual ~RecBase(){}

  /*! \brief Returns the grid associated with this reconstruction */
  Grid & getGrid() 
  {return _grid;}

  /*! \brief Compute the reconstruction coefficients in a range of cells.
   * 
   * Do here all possible background work that does not require the 
   * knowledge of the point where the reconstruction will be needed.
   */
  virtual void compute(const DofVec<M> * u, int cellStart, int cellEnd) =0;

  /*! \brief Compute the reconstruction coefficients in the entire grid.
   * 
   * Do here all possible background work that does not require the 
   * knowledge of the point where the reconstruction will be needed.
   */
  virtual void compute(const DofVec<M> * u)
  { compute(u, _grid.beginOneGhost(), _grid.endOneGhost() ); }

  /*! \brief Eval method: actually evaluates the reconstruction for a
   * given cell and at a given relative location in (-0.5,0.5)
   */
  virtual t_u eval(int k, tC relPos) =0;

  /*! \brief Eval method: actually evaluates the reconstruction for a
   * given cell and at its left boundary
   */
  virtual t_u evalLeft(int k)
  { return eval(k, -0.5); }

  /*! \brief Eval method: actually evaluates the reconstruction for a
   * given cell and at its right boundary
   */
  virtual t_u evalRight(int k)
  { return eval(k, 0.5); }

  /*! \brief Returns a pointer to the reconstruction coefficients for a given cell
   *
   * Warning: if the rec has no coefficients (e.g. p-wise constant), then
   * it returns 0. Just be careful when using the returned pointer...
   */
  virtual const t_u * getCoeffRec(int k) const =0;

  /*! \brief Override the reconstruction coefficients for a given cell
   *
   * Folks, don't try this at home!
   * ... and don't complain if you mess it up.
   *
   * This copies the reconstruction coefficients from the given pointer.
   * It assumes that the correct number of t_u's can be read at that pointer.
   */
  virtual void overrideCoeffRec(int k, const t_u* coeff)=0;

  /*! \brief Returns the number of ghosts needed by the reconstruction
   * procedure */
  virtual int needsGhosts() const =0;

  //! \brief Returns the theoretical order of accuracy of the recontruction procedure
  virtual int getOrder() const =0;

  /*! \brief Update pointer to data and resize storage for coefficients
   * 
   * The method compute ought to call this method first and then resize
   * any other additional storage needed during the computation, if any.
   */
  virtual void setUAvgAndResize(const DofVec<M> * u) =0;

  /*! Create another reconstruction of the same type
   *
   * Override this in derived classes with a method that creates a Rec object of the same
   * type of the derived class and returns a pointer to it.
   *
   * This is similar to a "clone" method, but we allow to attach the new reconstruction
   * to a different grid.
   */
  virtual RecBase<M> * duplicate(Grid & grid/*! grid to which the new rec will be attached */ ) const = 0;
  //! Prints info on the reconstruction
  virtual void print() const =0 ;

protected:
  Grid & _grid; //!< reference to the grid
};

#endif

/** @} */
