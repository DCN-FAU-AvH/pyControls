#pragma once

#include "config.h"
#include "utils/parseOptions.hpp"
#include "grid/grid.hpp"
#include "grid/gridutils.hpp"
#include "dof/dofvector.hpp"
#include "utils/save2file.hpp"
#include "utils/gaussRules.hpp"
#include "reconstructions/constRec.hpp"
#include "reconstructions/linRec.hpp"
#include "reconstructions/centralRec.hpp"
#include "reconstructions/cweno.hpp"
#include "claws/claw_psystem.hpp"
// templates

#include "claws/claw_burgers.tpp"
#include "claws/claw_linear.tpp"
#include "claws/claw_gas.tpp"
#include "claws/claw_swe.tpp"
#include "numfluxes/llfFlux.hpp"
#include "numfluxes/hllcFlux.hpp"
#include "timeintegration/timeEuler.hpp"
#include "timeintegration/erk.hpp"

struct cweno3_psys{
  explicitRungeKutta<3> * timeIntegrator_ptr;
  CLawPsysBCHandler * BC_ptr;
  CLawPsys clawObj;
  Grid * G_ptr;
  RecBase<3> * Rec_ptr;
  FluxBase<3> * numFlux_ptr;
  numSource<3> * Source_ptr;
  semidiscreteRHS<3> * doRHS;
  int cw_order = 3;
  cweno3_psys(double alpha, double gamma, double theta, double L, unsigned nr_cells,
		     unsigned  nr_ghosts, double epsilon,
	      unsigned alphaPower) : clawObj(CLawPsys(alpha, gamma)){
    G_ptr = new Grid(0.0, L, nr_cells, nr_ghosts);
    G_ptr->setLeftBC(Grid::BC_FREEFLOW);
    G_ptr->setRightBC(Grid::BC_FREEFLOW);
    
    Rec_ptr = new cweno<3>(*G_ptr, cw_order, cweno<3>::CW, epsilon , alphaPower);
    numFlux_ptr = new llfFlux<CLawPsys>(clawObj);
    Source_ptr = new PsysFrictionSource(theta, [] (double x) {return 0.0;}, *G_ptr, *Rec_ptr);
    BC_ptr = new CLawPsysBCHandler();

    doRHS = new naifRHS<3>(*G_ptr, *Rec_ptr, *numFlux_ptr, *BC_ptr , *Source_ptr);
    const int ghosts = doRHS->needsGhosts();
    G_ptr->setGhosts(ghosts);
    ExplicitButcherTableaux * erkTableaux = &getExplicitButcherTableaux<ERKtype>(ERK_ssp3);
    timeIntegrator_ptr = new explicitRungeKutta<3>(*erkTableaux, *G_ptr ,*doRHS);
    
  }
  ~cweno3_psys(){
    delete G_ptr;    
  delete Rec_ptr;
  delete numFlux_ptr;
  delete timeIntegrator_ptr;
  if (Source_ptr)
    delete Source_ptr;
  delete doRHS;
  }
  explicitRungeKutta<3> get_timeIntegrator(){
    return *timeIntegrator_ptr;
  }
  Grid get_Grid(){
    return *G_ptr;
  }

   void set_bc(bool left, std::vector<double> U_std){
      DoF<3, tD> U;
      for(auto j = 0; j < CLawPsys::M; ++j)
         U[j] = U_std[j];
      if (left)
        BC_ptr->left_bc =  std::vector<DoF<3, tD>>(cw_order, U);
      else
        BC_ptr->right_bc =  std::vector<DoF<3, tD>>(cw_order, U);
    }
    void set_both_bc(std::vector<double> left, std::vector<double> right){
      set_bc(true, left);
      set_bc(false, right);
    }
};

