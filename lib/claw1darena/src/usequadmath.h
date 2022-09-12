/*
  This file is part of claw1d, a software library to discretize
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

#ifndef USEQUADMATH_H
#define USEQUADMATH_H

#include <quadmath.h>
#include <iostream>

namespace std{
#ifndef STD_HAVE_ABSQUAD
  inline __float128 abs( __float128 x )
    { return fabsq( x ); }
#endif

  inline __float128 sqrt( __float128 x )
    { return sqrtq( x ); }

  inline __float128 sin( __float128 x )
    { return sinq( x ); }

  inline __float128 cos( __float128 x )
    { return cosq( x ); }

  inline __float128 exp( __float128 x )
    { return expq( x ); }

  inline __float128 pow( __float128 x , __float128 y )
    { return powq( x , y ); }

  inline __float128 max( __float128 x , __float128 y )
    { return fmaxq( x , y ); }

  inline std::ostream& operator<< (std::ostream& os, const __float128& f) {
    char* y = new char[100];
    quadmath_snprintf(y, 100, "%.7Qe", f) ;
    os.precision(30);
    os<<y;
    delete[] y;
    return os;
  }
}

#endif
