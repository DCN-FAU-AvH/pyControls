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

// class for efficient fixed size vectors,
// a sort of replacement of blitz::tinyvec
// This file contains the templated implementations.

/** @addtogroup dof
 *
 * @{
 */

/*! \file dof.tpp
 *  \brief Definitions of DoF<M> mambers
 */

#include "dof.hpp"
#include <cmath>

template<int N,class T>
inline DoF<N,T> operator*(const T &d, const DoF<N,T>& u)
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = d * u._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> operator/(const DoF<N,T>& u , const T &d)
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = u._data[i] / d;
  return r;
}

template<int N,class T>
inline DoF<N,T> operator/(const T &d, const DoF<N,T>& u)
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = d / u._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> operator-(const DoF<N,T>& u)
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = -u._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> operator-(const DoF<N,T>& u , const DoF<N,T>& v )
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = u._data[i] - v._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> operator+(const DoF<N,T>& u , const DoF<N,T>& v )
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = u._data[i] + v._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> operator+(const DoF<N,T>& u , const T& d )
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = u._data[i] + d;
  return r;
}

template<int N,class T>
inline DoF<N,T> eWiseProduct(const DoF<N,T>& u , const DoF<N,T>& v )
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = u._data[i] * v._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> eWiseDivide(const DoF<N,T>& u , const DoF<N,T>& v )
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = u._data[i] / v._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> eWiseMax(const DoF<N,T>& u , const DoF<N,T>& v )
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = std::max(u._data[i], v._data[i]);
  return r;
}

template<int N,class T>
inline DoF<N,T> eWiseAbs(const DoF<N,T>& u)
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = std::abs(u._data[i]);
  return r;
}

template<int N,class T>
inline DoF<N,T> eWiseSquare(const DoF<N,T>& u)
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = u._data[i] * u._data[i];
  return r;
}

template<int N,class T>
inline DoF<N,T> eWisePow(const DoF<N,T>& u, const tD& exp)
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = std::pow(u._data[i] , exp);
  return r;
}

template<int N,class T>
inline DoF<N,T> eWiseDiv(const T &d , const DoF<N,T>& v )
{
  DoF<N,T> r;
  for (uint i=0; i<N; i++)
    r._data[i] = d / v._data[i];
  return r;
}

template<int N,class T>
std::ostream& operator<<(std::ostream& s, const DoF<N,T> u) {
  s << "(";
  for (int i=0; i<N; ++i)
    s << u._data[i] << (i<N-1?", ":"");
  s <<  ")" ;
  return s ;
}

template<int N,class T>
std::istream& operator>>(std::istream& s, DoF<N,T> &u) {
  for (int i=0; i<N; ++i)
    s >> u._data[i];
  return s ;
}

template<int N, class T>
T max(DoF<N,T>& u){
  T M=u[0];
  for (uint i=1; i<N; ++i){
    if (M<u[i])
      M=u[i];
  }
  return(M);
}

template<int N, class T>
T min(DoF<N,T>& u){
  T m=u[0];
  for (uint i=1; i<N; ++i){
    if (m>u[i])
      m=u[i];
  }
  return(m);
}

/** @} */
