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

#ifndef __SRCCONFIG_HH
#define __SRCCONFIG_HH

#include "Claw1dArena_config.h"

//#ifdef DATATYPE_DOUBLE
typedef double tC; //type for coordinates
typedef double tD; //type for data (cell averages, ...)
//#define TD_IS_DOUBLE
//#endif

#ifdef DATATYPE_FLOAT128
#include "usequadmath.h"
typedef __float128 tC; //type for coordinates
typedef __float128 tD; //type for data (cell averages, ...)
#endif

#ifdef DATATYPE_CUSTOM
typedef float  tC; //type for coordinates
typedef double tD; //type for data (cell averages, ...)
#define TD_IS_DOUBLE
#endif

#endif
