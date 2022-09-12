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

/*! \file constRec.hpp
 *  \brief Declaration of ConstRec
 */

#ifndef CONSTREC_HH
#define CONSTREC_HH

#include "recBase.hpp"

/*! \brief Piece-wise constant recontruction
 *
 * This specialization of Recbase simply performs a piece-wise constant
 * reconstruction. Templated on the number of conserved variables.
 *
 * When compute is called, a pointer to the vector of the conserved
 * variables is stored internally and later used by the eavl methods.
 */
template <int M>
class ConstRec: public RecBase<M>{
public:
  typedef DoF<M, tD> t_u;

  ConstRec(Grid & grid);

  virtual void compute(const DofVec<M> * u, int, int);

  virtual t_u eval(int k, tC);

  virtual const t_u * getCoeffRec(int k) const;

  virtual void overrideCoeffRec(int k, const t_u* coeff);

  virtual int needsGhosts() const;

  virtual int getOrder() const;

  virtual void setUAvgAndResize(const DofVec<M> * u)
    { _uData=u; }

  virtual ConstRec<M> * duplicate(Grid & grid) const {
    std::cout<<"*** Called duplicate method of ConstRec ***\n";
    return new ConstRec<M>(grid);
  }

  //! Prints info on the reconstruction
  virtual void print() const {
    std::cout << "piece-wise constant reconstruction" << std::endl;
  };

private:
  DofVec<M> const * _uData;

};

#endif

/** @} */
