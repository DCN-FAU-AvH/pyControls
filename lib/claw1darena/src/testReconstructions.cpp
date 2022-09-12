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
#include "dof/dofvector.hpp"
#include "utils/gaussRules.hpp"
#include "utils/testrec.hpp"
#include "reconstructions/constRec.hpp"
#include "reconstructions/linRec.hpp"
#include "reconstructions/centralRec.hpp"
#include "reconstructions/cweno.hpp"

int main(int argc, char ** argv){
  std::cout << "===============================================" << std::endl;
  std::cout << "|            Test reconstructions             |" << std::endl;
  std::cout << "|                   "
            << "ver "<<claw1dArena_VERSION_MAJOR
            << "."   <<claw1dArena_VERSION_MINOR
            << "                   |" << std::endl;
  std::cout << "===============================================" << std::endl;

  Options opzioni;
  opzioni.parseCmdLine(argc, argv);

  const int N = opzioni.getNcells();
  const tC dx = 1.0/N;

  std::cout << "Using " << N << " cells (dx= " << dx << ")\n";
  std::cerr << dx << ", ";

  TestRec testObj;
  const TestRec::PBType pbType = testObj.getPBType(opzioni.getPBString());

  //Grid of 1 cell, with critical point in the middle
  const tC gridCenter = testObj.getCPLocation(pbType);
  Grid G (gridCenter-0.5*dx,gridCenter+0.5*dx,1);
  std::cout << "Critical point at " << gridCenter <<std::endl;

  //Choose reconstruction type from command line options
  RecBase<1> * Reconstruction;
  if (opzioni.getRecName().compare("const")==0){
    Reconstruction = new ConstRec<1>(G);
  }
  else if (opzioni.getRecName().compare("linear")==0){
    Reconstruction = new LinRec<1>(G);
  }
  else if (opzioni.getRecName().substr(0,7).compare("central")==0){
    std::string centralOrder(opzioni.getRecName().substr(7));
    try{
      int order = std::stoi(centralOrder);
      Reconstruction = new CentralRec<1>(G,order);
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
      cweno<1> * tmpRec;
      tmpRec = new cweno<1>(G,order,cweno<1>::CWZdb, epsilon, alphaPower);
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
          epsPower=3;
          break;
        default:
          epsPower=3;
        }
      }
      if (epsPower>0)
        epsilon = std::pow(G.getDx(),epsPower);
      else
        epsilon = std::pow(10.0 , epsPower);
      cweno<1> * tmpRec;
      tmpRec = new cweno<1>(G,order,cweno<1>::CWZ, epsilon, alphaPower);
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
      Reconstruction = new cweno<1>(G,order,cweno<1>::CW,epsilon,alphaPower);
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

  std::cout << "Reconstructions: "; Reconstruction->print();

  const int ghosts = Reconstruction->needsGhosts();
  G.setGhosts(ghosts);
  std::cout << "Setting up grid with " << ghosts << " ghosts" <<std::endl;

  DofVec<1> U;
  U.resize(G.getFullSize());

  //Quadrature rule for the initial data
  gaussRule quadrature;
  int initQuadOrder = opzioni.getQuadOrder();
  if (initQuadOrder<1)
    initQuadOrder = 1 + Reconstruction->getOrder();
  quadrature.chooseOrder(initQuadOrder);
  std::cout << "Quadrature rule for initial datum: Gauss-Legendre with "
    << quadrature.getNNodes() << " node(s)." << std::endl;

  std::cout << "Test: " << testObj.getPBName(pbType) <<std::endl;
  for (int i=G.beginFull(); i<G.endFull(); ++i){
    tC x = G.getXCenter(i);
    testObj.setU0(pbType, x, dx, quadrature, U[i]);
  }

  Reconstruction->compute(&U , G.beginPhysical(), G.endPhysical());

  const int Nrec=1;
  DoF<1,tD> err(0.);
  DoF<1,tD> errInf(0.);
  for (int i=-Nrec; i<=Nrec; i++){
    tC tRec = i/(2.0*Nrec);
    tC xRec = G.getXCenter(G.beginPhysical()) + tRec*dx;
    DoF<1,tD> uRec = Reconstruction->eval(G.beginPhysical(), tRec);
    DoF<1,tD> uExa;
    testObj.U0(pbType,xRec,uExa);
    err = eWiseAbs(uRec - uExa);
    errInf = eWiseMax(errInf , err);
  }
  std::cout << "Errore 0+1/2: " << err
    << " Errore in norma infinito " << errInf
    << std::endl;
  std::cerr << err[0] << ", "
    << errInf[0] << std::endl;
  delete Reconstruction;

  return 0;
}
