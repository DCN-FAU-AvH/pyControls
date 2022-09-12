#include "claw1darena_all.hpp"


void set_bc(CLawPsysBCHandler & BC, RecBase<3> * Rec_ptr, int cw_order, std::vector<double> left, std::vector<double> right){
  for(auto i = 0; i < Rec_ptr->needsGhosts(); ++i){
      BC.left_bc.resize(cw_order/2);
      BC.right_bc.resize(cw_order/2);
    for(auto j = 0; j < CLawPsys::M; ++j){
      BC.left_bc[i][j] = left[j];
      BC.right_bc[i][j] = right[j];
    }
  }
}
 

int main(){
  typedef CLawPsys tClaw;
  double  L = 100.0;
  int nr_cells = 100;
  double epsilon = 1.0e-20;
  int alphaPower = 0;
  int cw_order = 3;
  double _alpha = 1.0;
  double _gamma = 1.0;
  double _theta = 0.7;
  auto nr_ghosts = cw_order;


// CLawPsys clawObj(_alpha, _gamma);
// Grid G = Grid(0.0, L, nr_cells, nr_ghosts);
// // G.setLeftBC(Grid::BC_DIRICHLET);
// // G.setRightBC(Grid::BC_DIRICHLET);

// G.setLeftBC(Grid::BC_FREEFLOW);
// G.setRightBC(Grid::BC_FREEFLOW);

 
// RecBase<3> * Rec_ptr = new cweno<3>(G, cw_order, cweno<3>::CW, epsilon , alphaPower);
// FluxBase<3> * numFlux_ptr = new llfFlux<CLawPsys>(clawObj);
// numSource<3> * Source_ptr = new ZeroSource<3>(G);
// CLawPsysBCHandler BC;
// BCHandler<3> * BC_ptr = &BC;
// semidiscreteRHS<3> * doRHS = new naifRHS<3>(G, *Rec_ptr, *numFlux_ptr, *BC_ptr , *Source_ptr);
// const int ghosts = doRHS->needsGhosts();
// G.setGhosts(ghosts);
// ExplicitButcherTableaux * erkTableaux = &getExplicitButcherTableaux<ERKtype>(ERK_ssp3);
  // auto timeIntegrator = explicitRungeKutta<3>(*erkTableaux, G ,*doRHS);

  auto cw = cweno3_psys(_alpha, _gamma, _theta, L, nr_cells, nr_ghosts, epsilon, alphaPower);
 auto timeIntegrator_ptr  = cw.timeIntegrator_ptr;
 auto G_ptr = cw.G_ptr;
 timeIntegrator_ptr->setCFL(0.4);
 DofVec<tClaw::m> U(0), U1(0);
 U.resize(G_ptr->getFullSize());
 U1.resize(G_ptr->getFullSize());

 //set initial
 for(auto i = G_ptr->beginPhysical(); i < G_ptr->endPhysical(); ++i ){
   auto x = G_ptr->getXCenter(i);
    if( x < 50.0)
      U[i][0] = 1.0;
    else
      U[i][0] = 0.125;
    U[i][1] = 0.0;
    U[i][2] = 0.0;
 }

 double t = 0;
 auto Tfinal = 5.0;
 cw.set_both_bc({1.0, 0.0, 0.0}, {0.125, 0.0, 0.0});
while (t < Tfinal){
   timeIntegrator_ptr->setDtMax(Tfinal - t);
   std::cout << t << std::endl;
   auto dt = timeIntegrator_ptr->advance(U,t, U1);
   t = t + dt;
   std::swap(U, U1);
 }
 return 0;
}
