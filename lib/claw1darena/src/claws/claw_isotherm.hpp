#pragma once


#include "../grid/grid.hpp"
#include "clawbase.hpp"
#include "bcHandler.hpp"
#include "num_source.hpp"
#include "../reconstructions/recBase.hpp"
#include "../dof/dofvector.hpp"
#include "../utils/gaussRules.hpp"
#include "../config.h"

//! \brief CLaw class defining Euler gas dynamics equation
class CLawIsotherm : public CLawBase<CLawGas,3>{
public:
  enum {M=3}; //!< enum for number of conserved quantities
  typedef CLawBase<CLawGas,3> t_CLawBase;
  //!< type for a set of conserved quantities defining a gas state
  typedef typename t_CLawBase::t_u t_u;

  CLawIsotherm(tD c2);

    typedef DoF<M,tD> t_u;           //!< type for a set of conserved quantities
  typedef std::array<std::string,M> t_names;

  CLawIsotherm() {} //!< Default constructor
  ~CLawIsotherm() {} 

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
