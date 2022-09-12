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

/*! \file centralRec.hpp
 *  \brief Declaration of CentralRec
 */

#ifndef CENTRALREC_HH
#define CENTRALREC_HH

#include "recBase.hpp"
#include "udivdiff.hpp"

/*! \brief A central unlimited reconstruction of arbitary (odd) order
 *
 * This reconstruction makes use of the CENTRAL UNLIMITED polinomial.
 * It will thus be UNSTABLE and is of some use only for checking the
 * order of methods on smooth problems...
 *
 * Templated on the number of conserved variables to be reconstructed.
 */
template <int M>
class CentralRec: public RecBase<M>{
public:
  typedef DoF<M, tD> t_u; //!< type for a set of conserved quantities

  CentralRec(Grid & grid, int accuracy=3);

  void setAccuracy(int accuracy);
  int  getAccuracy();

  virtual void compute(const DofVec<M> * u, int cellStart, int cellEnd);

  virtual t_u eval(int k, tC);

  virtual const t_u * getCoeffRec(int k) const;

  virtual void overrideCoeffRec(int k, const t_u* coeff);

  virtual int needsGhosts() const;

  virtual int getOrder() const;

  virtual void setUAvgAndResize(const DofVec<M> * u);

  virtual CentralRec<M> * duplicate(Grid & grid) const {
    std::cout<<"*** Called duplicate method of CentralRec ***\n";
    return new CentralRec<M>(grid);
  }

  //! Prints info on the reconstruction
  virtual void print() const {
    std::cout << "UNLIMITED central reconstruction of order " << _accuracy << std::endl;
  };

private:
  int _accuracy;            //!< degree of the central polynomial
  UDivDiff<M> _udd;         //!< UnDividedDifferences object for polynomial computations
  std::vector<std::vector<t_u>> _poly;     //!< storage for the reconstruction polynomials
};

#endif

/** @} */
