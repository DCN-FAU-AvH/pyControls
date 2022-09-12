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

#include "polyTools.hpp"
#include <math.h>

/** @addtogroup reconstructions
 *
 * @{
 */

/*! \file polyTools.tpp
 *  \brief Implementation of polyTools functions
 */


/*! \brief Compute the Jiang-Shu indicator
 *
 * The coefficients of the polynomial should be in the basis ((x-x_j)/h)^k,
 * stored starting from the least significant one, like in the output
 * of UDiffDiv::interPoly().
 */
template <int M>
DoF<M,tD> computeIndicator(int deg/*!degree*/, DoF<M,tD> *P/*!polynomial coefficients*/){
  tD quadBase;
  DoF<M,tD>  Q[2*deg-1];
  for(int k=0;k<2*deg-1;++k)
    Q[k]=DoF<M,tD>(0.0);
  DoF<M,tD>  p[deg+1];
  for(int k=0;k<deg+1;++k)
    p[k]=DoF<M,tD>(0.0);
  
  DoF<M,tD> I(0.0);
  
  /* Change the ascending order of the coefficients */
  for (int j=0;j<deg+1;++j)
    p[j] = P[deg-j];
    /* Loop on the derivative */
    for (int d=deg;d>0;--d){
      /* Differentiate */
      for (int k=deg;k>0;--k)
        p[k] = tD((deg-k+1))*p[k-1];
        /* Set Q to 0 */
        for (int j=0;j<2*d-1;++j)
          Q[j]=DoF<M,tD>(0.0);
          /* Compute the square of the polynomial */
          for (int j=0;j<d;++j)
            for (int i=0;i<d;++i)
              Q[2*d-2-i-j] += eWiseProduct(p[deg-i],p[deg-j]);
              /* Update I */
              for (int j=0;j<d;++j){
                quadBase = pow(2.0,-2*j) / (2*j+1);
                I += quadBase * Q[2*(d-j)-2];
              }
    }

  return I;
}

/*! \brief Evaluate a polynomial P
 *
 * The coefficients should be with respect to the basis (x-x_j)^k/h
 *
 * The value pos should thus be between -0.5 and 0.5.
 */
template <int M>
DoF<M,tD> evalPoly(int deg/*!degree*/, tC pos/*!evaluation point*/, DoF<M,tD> *P/*!polynomial coefficients*/){
  DoF<M,tD> value=DoF<M,tD>(0.0);
  
  value = tD(pos)*P[deg];
  
  for(int k=deg-1;k>0;--k)
    value = tD(pos)*(P[k]+value);
    
  value = P[0]+value;
    
  return value;
}

/** @} */
