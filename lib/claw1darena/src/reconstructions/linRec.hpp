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

/** @addtogroup reconstructions
 *
 * @{
 */

/*! \file linRec.hpp
 *  \brief Declaration of LinRec
 */

#ifndef LINREC_HH
#define LINREC_HH

#include "recBase.hpp"
#include "udivdiff.hpp"

/*! \brief Implements slope-limited linear reconstruction
 *
 * Implements a minmod reconstruction.
 *
 * Templated on the number of conserved variables to be reconstructed.
 */
template <int M>
class LinRec: public RecBase<M>{
public:
  typedef DoF<M, tD> t_u; //!< type for a set of conserved quantities

  LinRec(Grid & grid);

  virtual void compute(const DofVec<M> * u, int cellStart, int cellEnd);

  virtual t_u eval(int k, tC);

  virtual const t_u * getCoeffRec(int k) const;

  virtual void overrideCoeffRec(int k, const t_u* coeff);

  virtual int needsGhosts() const;

  virtual int getOrder() const;

  virtual void setUAvgAndResize(const DofVec<M> * u);

  virtual LinRec<M> * duplicate(Grid & grid) const {
    std::cout<<"*** Called duplicate method of LinRec ***\n";
    return new LinRec<M>(grid);
  }

  //! Prints info on the reconstruction
  virtual void print() const {
    std::cout << "piece-wise linear reconstruction, using minmod" << std::endl;
  };

private:
  DofVec<M> const * _uData; //!< pointer to the cell averages of the data
  DofVec<M> _slopeData;     //!< storage for the slopes of the reconstruction
  DofVec<M> _uSlope;        //!< object to store the unlimited slopes

};

#endif

/** @} */
