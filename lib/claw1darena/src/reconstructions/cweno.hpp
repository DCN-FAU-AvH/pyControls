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

/*! \file cweno.hpp
 *  \brief Declaration of CWENO
 */

#ifndef CWENO_HH
#define CWENO_HH

#include "recBase.hpp"
#include "udivdiff.hpp"

/*! \brief Implements central-weno reconstruction
 * 
 * Templated on the number of conserved variables to be reconstructed.
 */
template <int M>
class cweno: public RecBase<M>{
public:
  typedef DoF<M, tD> t_u; //!< type for a set of conserved quantities

  enum CWType{CW/*!CWENO of Math of Comp paper*/,
             CWZdb/*!CWENOZ of Comp. and Fluids paper, using Don-Borges optimal tau*/,
             CWZ/*!The new CWENOZ*/};

  cweno(Grid & grid,int accuracy=3,CWType type=CW,tD epsilon=1.0e-12, int alphaPower=0);

  void setAccuracy(int accuracy);
  int getAccuracy();
  
  void setType(CWType type);
  CWType getType() const;
  
  void setAlphaPower(int alphaPower);
  int  getAlphaPower() const;
  void setDefaultAlphaPower();

  void setEpsilon(tD epsilon);
  tD getEpsilon() const;

  void setIP0(const bool IP0=false);
  bool getIP0() const;
  bool getDefaultIP0() const;

  virtual void compute(const DofVec<M> * u, int cellStart, int cellEnd);

  virtual t_u eval(int k, tC);

  virtual const t_u * getCoeffRec(int k) const;

  virtual void overrideCoeffRec(int k, const t_u* coeff);

  virtual int needsGhosts() const;

  virtual int getOrder() const;

  virtual void setUAvgAndResize(const DofVec<M> * u);

  virtual cweno<M> * duplicate(Grid & grid) const {
    //std::cout<<"*** Called duplicate method of cweno ***\n";
    cweno<M> * tmpRec = new cweno<M>(grid, _accuracy, _type, _epsilon);
    tmpRec->setIP0(_useIP0);
    return tmpRec;
  }

  //! Prints info on the reconstruction
  virtual void print() const;

private:
  CWType _type;
  int _accuracy; //!< accuracy of the reconstruction
  tD _epsilon; //!< epsilon for the computation of nonlinear weights
  int _alphaPower; //!< power used in the computation of nonlinear weights
  t_u _tau; //!< tau for the computations of nonlinear weights in CWZ
  int _nPolyTot; //!< total number of polynomials
  int _nPolyLow; //!< number of low degree polynomials
  int _degLow; //!< degree of low degree polynomials
  bool _useIP0; //!< choice of I0 in CW and CWZ: if true, use I0=I[P0], else I0=I[Popt]

  UDivDiff<M> _udd; //!< UnDividedDifferences object for polynomial computations

  std::vector<tD> _optWeights; //!< vector of optimal weights
  std::vector<t_u> _w; //!< vector of nonlinear weights
  std::vector<std::vector<t_u>> _polyHigh; //!< storage for the reconstruction polynomials
  std::vector<std::vector<t_u>> _polyLow; //!< storage for the low degree polynomials
  std::vector<t_u> _polyInds; //!< storage for indicators

  void setIdealWeights(int accuracy);
};

#endif

/** @} */
