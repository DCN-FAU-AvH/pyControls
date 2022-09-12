#pragma once
#include <cassert>
#include "../grid/grid.hpp"
#include "../grid/gridutils.hpp"
#include "clawbase.hpp"
#include <functional>
#include "bcHandler.hpp"
#include "num_source.hpp"

//class CLawPsysBCHandler;
class CLawPsys : public CLawBase<CLawPsys,3>{
public:
  enum {M=3}; //!< enum for number of conserved quantities
  enum {RHO=0, MOM=1};
  enum PBType {DUMMY=0};
  typedef CLawBase<CLawPsys,3> t_CLawBase;
  //!< type for a set of conserved quantities defining a gas state
  typedef typename t_CLawBase::t_u t_u;
  tD alpha, gamma;
  CLawPsys(tD _alpha, tD _gamma) : alpha(_alpha), gamma(_gamma) {};


  t_u getFlMax(const t_u & u, tC, tC, tD &lMax) const{
    assert( u[RHO] >=0 );
    const tD v = u[MOM]/u[RHO];
    const tD p = alpha*pow(u[RHO], gamma);
    assert ( p >=0 );
    const tD c= std::sqrt( gamma * p / u[RHO]);
    lMax = std::abs(v)+c;
    t_u f;
    f[RHO] = u[MOM];
    f[MOM] = v * u[MOM] + p;
    f[2]  = 0.0; //todo delete when refactoring to M=2 
    return f;
  }

  tD getLambdaMax(const t_u & u, tC , tC) const{
    assert( u[RHO] >=0 );
    const tD v = u[MOM]/u[RHO];
    const tD p =  alpha*pow(u[RHO], gamma);
    assert ( p >=0 );
    const tD c= std::sqrt( gamma * p / u[RHO]);
    return (std::abs(v)+c);
  }

  t_u getF(const t_u & u, tC, tC) const{
    assert( u[RHO] >=0 );
    const tD v = u[MOM]/u[RHO];
    const tD p = alpha*pow(u[RHO], gamma);
    assert ( p >=0 );
    const tD c= std::sqrt( gamma * p / u[RHO]);
    t_u f;
    f[RHO] = u[MOM];
    f[MOM] = v * u[MOM] + p;
    f[2]  = 0.0; //todo delete when refactoring to M=2 
    return f;
  }

    //! \brief Converts a string into a PBType
  int getPBType(const std::string& PBString) const{
    throw std::runtime_error("getPBtype not implemented!");
  }

  //! \brief Converts a PBType into a string
  std::string getPBName(int PBType) const{
    throw std::runtime_error("getPBName not implemented!");
  }
  
//! \brief Returns the function that represents the boundary data
  template <typename T>
void setBCFuncs(CLawPsys::PBType pb,T & bcH){
  // bcH.bcFunLeft  = NULL;
  // bcH.bcFunRight  = NULL;
}


};

// todo add Romberg source like in SWE
class PsysFrictionSource : public numSource<3>{
public:
  typedef numSource<3>::t_u t_u; //!< type for a set of conserved quantities

  //! Constructor from grid and reconstruction object
  PsysFrictionSource(double _theta,  std::function<double(double)> _w,
		     Grid & grid, RecBase<3> & rec/*! the reconstruction object*/):
    numSource<3>(grid),
    _rec(rec),
    theta(_theta),
    w(_w)
  {}

  ~PsysFrictionSource(){}

   t_u getSource(int idx, tC time) const{
    t_u S(0.);
    double dx = _grid.getDx();
    double x = idx*dx; // todo for w(x)
    t_u recC = _rec.eval(idx  ,0.);
    t_u recL = _rec.eval(idx-1,0.);
    t_u recR = _rec.eval(idx+1,0.);

    S[1] = -0.5 * (theta + w(x))* recC[1]*fabs(recC[1])/recC[0];
    return S;
  }

  void setOrder(int order) {};
  int getOrder() const {return 1;};

private:
  double theta;
  std::function<double(double)> w;
  //  gaussRule qRule; //!< reference to the object of gaussRule
  RecBase<3> & _rec; //!< reference to the object RecBase
};


class CLawPsysBCHandler : public BCHandler<3>{
public:
  CLawPsysBCHandler(){
    DoF<3, tD> left_U, right_U;
    left_U[0] = 1.0;
    left_U[1] = 0.0;
    left_U[2] = 0.0;
    right_U[0] = 0.125;
    right_U[1] = 0.0;
    right_U[2] = 0.0;
    left_bc =  std::vector<DoF<3, tD>>(3, left_U);
    right_bc =  std::vector<DoF<3, tD>>(3, right_U);    
  }
  void setLeftGhosts  (Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U)  {
    if (bc == Grid::BC_PERIODIC)
       setLeftPeriodicGhosts<3>(grid , U);
    // else if (bc == Grid::BC_DIRICHLET){ //maybe later
    //   DoF<M,tD> v = bc_left[0];
    //   setLeftDirichletGhosts<3>(grid, U, v);
    // }
    else {
      int k = 0 ;
      for (int i=grid.rendLeftGhosts(); i<=grid.rbeginLeftGhosts() ; i++){
        U[i] = left_bc[k];
        k++;
      }
    }
  }
  void setRightGhosts (Grid &grid, Grid::BCType bc, tC t, DofVec<3> &U){
    if (bc == Grid::BC_PERIODIC)
      setRightPeriodicGhosts<3>(grid , U);
    // else if (bc == Grid::BC_DIRICHLET){ //maybe later
    //   DoF<M,tD> v = bc_right[0];
    //   setRightDirichletGhosts<3>(grid , U, v);
    // }
    else {
      int k = 0;
      for (int i=grid.beginRightGhosts(); i<grid.endRightGhosts() ; i++){
        U[i] = right_bc[k];
        k++;
      }
    }
  }
  std::vector<DoF<3, tD>> left_bc;
  std::vector<DoF<3, tD>> right_bc;
private:
  bcFuncType<3> *bcFunLeft;
  bcFuncType<3> *bcFunRight;

  //friend void CLawPsys::setBCFuncs(PBType pb , CLawPsysBCHandler & bcH);
};

