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

/** @defgroup claws Representation of conservation laws
 *
 * @{
 */

/*! \file clawbase.hpp
 *  \brief Declaration of CLawBase
 */

#ifndef _CLAWBASE_HH
#define _CLAWBASE_HH

#include "../dof/dof.hpp"
#include "../config.h"
#include <string>

/*! \brief Base class to represent a system of M conservation law in 1d and its associated problems.
 *
 * The base class is templated on a class CLaw and on an integer M.
 * CRTP is employed and all methods are forwarded to the real implementation
 * which is contained in the template CLaw.
 * M is the number of conserved quantities.
 */
template<class CLaw, int M>
class CLawBase{
public:
  enum {m=M};                      //!< enum for number of conserved quantities
  typedef DoF<M,tD> t_u;           //!< type for a set of conserved quantities
  typedef std::array<std::string,M> t_names;

  CLawBase() {} //!< Default constructor
  virtual ~CLawBase() {} //!< virtual destructor

  //! \brief Returns the flux and the spectral radius on u, at x and t
  t_u getFlMax(const t_u& u, tC x, tC t, tD & lMax) const{
    return ReferToDerived().getFlMax(u,x,t,lMax);
  }

  //! \brief Returns the flux function on u, at x and t
  t_u getF(const t_u& u, tC x, tC t) const{
    return ReferToDerived().getF(u,x,t);
  }

  //! \brief Returns the spectral radius on u, at x and t
  tD getLambdaMax(const t_u& u, tC x, tC t) const{
    return ReferToDerived().getLambdaMax(u,x,t);
  }

  //! \brief Converts a string into a PBType
  int getPBType(const std::string& PBString) const{
    return ReferToDerived().getPBType(PBString);
  }

  //! \brief Converts a PBType into a string
  std::string getPBName(int PBType) const{
    return ReferToDerived().getPBName(PBType);
  }

  //! \brief Returns the names of the conserved quantities
  static const t_names & getVarnames(){
    return varnames;
  }
protected:
  const static std::array<std::string,M> varnames;  //!< names of the conserved quantities
private:
  //! \brief The CRTP trick
  const CLaw& ReferToDerived() const{
    return static_cast<const CLaw&>(*this);
  }
};

//! \brief Exception thrown when the name of a problem is not recognized
class UnknownProblemString{
public:
  UnknownProblemString(std::string name): _name(name) {}
  std::string _name; //!< offending string
};

//! \brief Exception thrown when a PBType is not recognized
class UnknownProblemId{
public:
  UnknownProblemId(int id): _id(id) {}
  int _id; //!< offending PBType
};

#endif

/** @} */
