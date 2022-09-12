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
/*
 * Note: all claw*.tpp need to be included even if only one of them is
 * actually needed in the final executable, otherwise a link error
 * occours due to the code in llfFlux.cpp and hllcFlux.cpp
 * See also the note in claw_gas.tpp for characteristic reconstructions
 * TODO!
*/
#include "claw_linear.tpp"
#include "claw_burgers.tpp"
#include "claw_gas.tpp"
#include "claw_swe.tpp"
#include "numfluxes/llfFlux.hpp"
#include "numfluxes/hllcFlux.hpp"
#include "timeintegration/timeEuler.hpp"
#include "timeintegration/erk.hpp"

int main(int argc, char ** argv){
  std::cout << "===============================================" << std::endl;
  std::cout << "|               Claw 1d Arena                 |" << std::endl;
  std::cout << "|                   "
            << "ver "<<claw1dArena_VERSION_MAJOR
            << "."   <<claw1dArena_VERSION_MINOR
            << "                   |" << std::endl;
  std::cout << "===============================================" << std::endl;

  Options opzioni;
  opzioni.parseCmdLine(argc, argv);

  std::cout << "Using " << opzioni.getNcells() << " cells\n";

#ifdef LINTRA
  typedef CLawLinear tClaw;
  typedef CLawLinearBCHandler tClawBCHandler;
#endif
#ifdef BURGERS
  typedef CLawBurgers tClaw;
  typedef CLawBurgersBCHandler tClawBCHandler;
#endif
#ifdef GASDIN
  typedef CLawGas tClaw;
  typedef CLawGasBCHandler tClawBCHandler;
#endif
#ifdef SWE
  typedef CLawSWE tClaw;
  typedef CLawSWEBCHandler tClawBCHandler;
#endif

  int N = opzioni.getNcells();
  tD cfl= opzioni.getCFL();

  tClaw clawObj;
  const tClaw::PBType pbType = clawObj.getPBType(opzioni.getPBString());
  tClawBCHandler clawBCHandlerObj;
  clawObj.setBCFuncs(pbType, clawBCHandlerObj);

  Grid G (clawObj.getLeftGridLimit (pbType),
          clawObj.getRightGridLimit(pbType),
          N);
  G.setLeftBC (clawObj.getBCLeft (pbType));
  G.setRightBC(clawObj.getBCRight(pbType));
  const tC dx = G.getDx();

  RecBase<tClaw::m> * Reconstruction;
  if (opzioni.getRecName().compare("const")==0){
    Reconstruction = new ConstRec<tClaw::m>(G);
  }
  else if (opzioni.getRecName().compare("linear")==0){
    Reconstruction = new LinRec<tClaw::m>(G);
  }
  else if (opzioni.getRecName().substr(0,7).compare("central")==0){
    std::string centralOrder(opzioni.getRecName().substr(7));
    try{
      int order = std::stoi(centralOrder);
      Reconstruction = new CentralRec<tClaw::m>(G,order);
    }
    catch (std::invalid_argument){
      std::cout << "Unknown reconstruction type " << opzioni.getRecName()
                << ". Aborting" <<std::endl;
      abort();
    }
  }
  else if (opzioni.getRecName().substr(0,8).compare("cwenozdb")==0){
    std::string cwenoOrder(opzioni.getRecName().substr(8));
    try{
      int order = std::stoi(cwenoOrder);
      int alphaPower = opzioni.getRecAlphaPower();
      int epsPower = opzioni.getRecEpsPower();
      tD epsilon =1.0e-14;
      if (epsPower==0){
        switch(order){
        case 3:
          epsPower=2;
          break;
        case 5:
          epsPower=4;
          break;
        case 7:
          epsPower=5;
          break;
        case 9:
          epsPower=5;
          break;
        default:
          epsPower=6;
        }
      }
      if (epsPower>0)
        epsilon = std::pow(G.getDx(),epsPower);
      else
        epsilon = std::pow(10.0 , epsPower);
      cweno<tClaw::m> * tmpRec;
      tmpRec = new cweno<tClaw::m>(G,order,cweno<tClaw::m>::CWZdb, epsilon, alphaPower);
      if (opzioni.forcedI0())
        tmpRec->setIP0(opzioni.getuseIP0());
      Reconstruction = tmpRec;
    }
    catch (std::invalid_argument){
      std::cout << "Unknown reconstruction type " << opzioni.getRecName()
                << ". Aborting" <<std::endl;
      abort();
    }
  }
  else if (opzioni.getRecName().substr(0,6).compare("cwenoz")==0){
    std::string cwenoOrder(opzioni.getRecName().substr(6));
    try{
      int order = std::stoi(cwenoOrder);
      int alphaPower = opzioni.getRecAlphaPower();
      int epsPower = opzioni.getRecEpsPower();
      tD epsilon =1.0e-14;
      if (epsPower==0){
        switch(order){
        case 3:
          epsPower=2;
          break;
        case 5:
          epsPower=4;
          break;
        default:
          epsPower=2;
        }
      }
      if (epsPower>0)
        epsilon = std::pow(G.getDx(),epsPower);
      else
        epsilon = std::pow(10.0 , epsPower);
      cweno<tClaw::m> * tmpRec;
      tmpRec = new cweno<tClaw::m>(G,order,cweno<tClaw::m>::CWZ, epsilon, alphaPower);
      if (opzioni.forcedI0())
        tmpRec->setIP0(opzioni.getuseIP0());
      Reconstruction = tmpRec;
    }
    catch (std::invalid_argument){
      std::cout << "Unknown reconstruction type " << opzioni.getRecName()
                << ". Aborting" <<std::endl;
      abort();
    }
  }
  else if (opzioni.getRecName().substr(0,5).compare("cweno")==0){
    std::string cwenoOrder(opzioni.getRecName().substr(5));
    try{
      int order = std::stoi(cwenoOrder);
      int alphaPower = opzioni.getRecAlphaPower();
      int epsPower = opzioni.getRecEpsPower();
      tD epsilon;
      if (epsPower==0)
          epsPower=2;
      if (epsPower>0)
        epsilon = std::pow(G.getDx(),epsPower);
      else
        epsilon = std::pow(10.0 , epsPower);
      Reconstruction = new cweno<tClaw::m>(G,order,cweno<tClaw::m>::CW,epsilon,alphaPower);
    }
    catch (std::invalid_argument){
      std::cout << "Unknown reconstruction type " << opzioni.getRecName()
                << ". Aborting" <<std::endl;
      abort();
    }
  }
  else{
    std::cout << opzioni.getRecName() << ": unknown reconstruction type" << std::endl;
    abort();
  }

  const bool useCharProj = opzioni.getCharProj();
  RecBase<tClaw::m> * cwiseReconstruction = NULL;
#if defined(GASDIN) && defined(TD_IS_DOUBLE)
  if (useCharProj){
    cwiseReconstruction = Reconstruction;
    Reconstruction = new GasCharRec(cwiseReconstruction);
  }
#else
  if (useCharProj){
    std::cout << "Sorry, I cannot do characteristic projection for this conservation law!" <<std::endl;
    abort();
  }
#endif

  FluxBase <tClaw::m> * numFlux;
  if (opzioni.getFluxName().compare("llf")==0){
    std::cout << "Numerical Flux: Local Lax-Friedrichs (Rusanov)" <<std::endl;
    numFlux = new llfFlux<tClaw> (clawObj);
  }
#ifndef SWE
  else if (opzioni.getFluxName().compare("hllc")==0){
    std::cout << "Numerical Flux: HLLC flux (Einfeldt's speed estimates)" <<std::endl;
    numFlux = new hllcFlux<tClaw> (clawObj);
  }
#endif
  else{
    std::cout << opzioni.getFluxName() << ": unknown flux type" << std::endl;
    abort();
  }

  typedef numSource<tClaw::M> tSource;
  tSource * Source=NULL;
#ifdef GASDIN
  switch (pbType){
  case CLawGas::PB_SOD_RADIAL_2D:
    Source = new CLawGasRadialSource2d(G, *Reconstruction);
    break;
  default:
    Source = new ZeroSource<3>(G);
  }
#endif
#ifdef LINTRA
    Source = new ZeroSource<1>(G);
#endif
#ifdef BURGERS
    Source = new ZeroSource<1>(G);
#endif

  semidiscreteRHS<tClaw::m> * doRHS = NULL;
#ifdef SWE
  const bool useWBhydro = opzioni.getWBhydro();
  if (useWBhydro){
    cwiseReconstruction = Reconstruction;
    std::cout << "Using hydrostatic well-balanced scheme" << std::endl;
    Reconstruction = new SWEhydroRec(Reconstruction);
    doRHS = new wbHydroRHS(G,*Reconstruction,*numFlux,clawBCHandlerObj);
  }
  else{
    //define source here!
    std::cout << "Naif unbalanced scheme" << std::endl;
    if ( Reconstruction->getOrder()==1)
      Source = new SWESourceUnb1(G, *Reconstruction);
    else
      Source = new SWESourceRomberg(G, *Reconstruction);
    doRHS = new naifRHS<3>(G,*Reconstruction,*numFlux,clawBCHandlerObj,*Source);
  }
#else
  doRHS = new naifRHS<tClaw::m>(G,*Reconstruction,*numFlux,clawBCHandlerObj,*Source);
#endif

  std::cout << "Reconstructions: "; Reconstruction->print();

  timeIntegration<tClaw::m> * timeIntegrator;
  if (opzioni.getTimeStepperName().compare("euler")==0){
    timeIntegrator= new Euler<tClaw::m>(G, *doRHS);
    std::cout << "Timestepper: Explicit Euler" <<std::endl;
  }
  else if (opzioni.getTimeStepperName().substr(0,3).compare("erk")==0)
  {
    std::string erkName(opzioni.getTimeStepperName().substr(3));
    ExplicitButcherTableaux * erkTableaux;
    try{
      int order = std::stoi(erkName);
      erkTableaux = &getExplicitButcherTableaux<int>(order);
    }
    catch (std::invalid_argument){
      if (erkName.compare("euler")==0)
        erkTableaux = &getExplicitButcherTableaux<ERKtype>(ERK_euler);
      else if (erkName.compare("heun")==0)
        erkTableaux = &getExplicitButcherTableaux<ERKtype>(ERK_heun);
      else if (erkName.compare("ssp3")==0)
        erkTableaux = &getExplicitButcherTableaux<ERKtype>(ERK_ssp3);
      else{
        std::cout << "Unknown ERK type " << erkName << ". Aborting" <<std::endl;
        abort();
      }
    }
    erkTableaux->print();
    timeIntegrator= new explicitRungeKutta<tClaw::m>(*erkTableaux,G, *doRHS);
  }
  else{
    std::cout << "Unknown timestepper method " << opzioni.getTimeStepperName()
      << ". Aborting" <<std::endl;
    abort();
  }

  timeIntegrator->setCFL(cfl);
  std::cout << " cfl=" << cfl << std::endl;

  const int ghosts = doRHS->needsGhosts();
  G.setGhosts(ghosts);
  std::cout << "Setting up grid with " << ghosts << " ghosts" <<std::endl;

  DofVec<tClaw::m> U,U1;
  U.resize(G.getFullSize());

  //Quadrature rule for the initial data
  gaussRule quadrature;
  int initQuadOrder = opzioni.getQuadOrder();
  if (initQuadOrder<1)
    initQuadOrder = std::max(Reconstruction->getOrder(), timeIntegrator->getOrder());
  quadrature.chooseOrder(initQuadOrder);
  std::cout << "Quadrature rule for initial datum: Gauss-Legendre with "
    << quadrature.getNNodes() << " node(s)." << std::endl;

  std::cout << "Test: " << clawObj.getPBName(pbType) <<std::endl;
  for (int i=G.beginPhysical(); i<G.endPhysical(); ++i){
    tC x = G.getXCenter(i);
    clawObj.setU0(pbType, x, dx, quadrature, U[i]);
  }
  save2File(G,U, (const char *) "initial.csv");

  tC tFinal(clawObj.getFinalTime(pbType));
  opzioni.getFinalTime(tFinal); //allow override from command line
  std::cout << " tFinal=" << tFinal <<std::endl;
  tC t(0.);
  tC nextOutput = tFinal/20.;

  int nStep(0);

  while (t<tFinal){
    timeIntegrator->setDtMax(tFinal-t);
    tC dt = timeIntegrator->advance(U,t,U1);
    t += dt;
    std::swap(U,U1);

    if (t>=nextOutput){
      std::cout << "Time " << t <<std::endl;
      nextOutput += tFinal/20.;
    }
    nStep++;
    //std::stringstream fileName;
    //fileName << "solution_" << nStep << ".csv";
    //save2File(G,U, fileName.str().c_str() );
  }

  std::cout<<"Computed " << nStep << " steps"<<std::endl;
  save2File(G,U, (const char *) "final.csv");

  //compute and print the error, if possible
  if (clawObj.haveExactSol(pbType)){
    for (int i=G.beginPhysical(); i<G.endPhysical(); ++i){
      tC x = G.getXCenter(i);
      clawObj.setUFinal(pbType, x, dx, quadrature, U1[i]);
    }
    tClaw::t_u normInf(0.), norm1(0.);
    tClaw::t_u TVERR(0.); //total variation
    for (int i=G.beginPhysical(); i<G.endPhysical(); ++i){
      tClaw::t_u diff = eWiseAbs(U[i]-U1[i]);
      normInf = eWiseMax( normInf , diff );
      norm1  += diff;
      if (i>G.beginPhysical())
        TVERR += eWiseAbs(U[i]-U[i-1]) - eWiseAbs(U1[i]-U1[i-1]);
    }
    norm1 *= dx;
    std::cout << "Inf-norm error: " << std::scientific << normInf << std::endl;
    std::cout << "  1-norm error: " << std::scientific << norm1 << std::endl;
    std::cerr << std::scientific << dx;
    for (int d=0; d<tClaw::m; ++d)
      std::cerr << " " << norm1[d];
    for (int d=0; d<tClaw::m; ++d)
      std::cerr << " " << normInf[d];
    for (int d=0; d<tClaw::m; ++d)
      std::cerr << " " << TVERR[d];
    std::cerr << std::endl;
  }

  delete Reconstruction;
  delete numFlux;
  delete timeIntegrator;
  if (Source)
    delete Source;
  if (cwiseReconstruction)
    delete cwiseReconstruction;
  delete doRHS;

  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "| If you use this software in your research,  |" << std::endl;
  std::cout << "| please cite it as                           |" << std::endl;
  std::cout << "| claw1dArena, by M. Semplice and G. Visconti |" << std::endl;
  std::cout << "| http://doi.org/10.5281/zenodo.2641725       |" << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

  return 0;
}
