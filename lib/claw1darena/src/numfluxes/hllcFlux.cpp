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

/** @addtogroup numfluxes
 *
 * @{
 */

/*! \file hllcFlux.hpp
 *  \brief Declaration of hllcFlux
 */

#include "hllcFlux.hpp"

//! \brief Constructor
hllcFlux<CLawGas>::hllcFlux(CLawGas & cLaw):
  _cLaw(cLaw)
  {}

/*! \brief The HLLC flux for the Euler gas dynamics equations
 *
 * The flux is computed following Toro's book, Section 10.6, but using
 * Einfeldt's speed estimates.
 * Note: the code for other estimates is currently in the source file
 * but commented out.
 */
typename hllcFlux<CLawGas>::t_u hllcFlux<CLawGas>::getF(const t_u& uL/*!left state*/, const t_u& uR/*!right state*/, tC, tC, tD & lL/*!set to speed of left-going wave in output*/, tD& lR/*!set to speed of left-going wave in output*/)
{
  t_u F;

  /* ==================================================
   * Implementation following Toro's book, Section 10.6
   * but using Einfeldt's speed estimates (see 10.5).
   * ==================================================
  */

  tD rhoL, vL, pL, aL;
  _cLaw.consToPrim(uL, rhoL, vL, pL, aL);
  tD rhoR, vR, pR, aR;
  _cLaw.consToPrim(uR, rhoR, vR, pR, aR);

  /* ==========================================
   *  Speed estimates: Toro's book Section 10.5
   * ==========================================
  */

  //Einfeldt's speed estimates: (10.52) -- (10.54)
  const tD sqrtRhoL = std::sqrt(rhoL);
  const tD sqrtRhoR = std::sqrt(rhoR);
  const tD eta2=0.5*sqrtRhoL*sqrtRhoR/((sqrtRhoL+sqrtRhoR)*(sqrtRhoL+sqrtRhoR));
  const tD d = std::sqrt((sqrtRhoL*aL*aL+sqrtRhoR*aR*aR)/
               (sqrtRhoL+sqrtRhoR) + eta2*(vR-vL)*(vR-vL));
  const tD sL = 0.5*(vL+vR)-d;
  const tD sR = 0.5*(vL+vR)+d;

  /*
  // Estimates using the Roe average: (10.49) -- (10.51)
  const tD sqrtRhoL = std::sqrt(rhoL);
  const tD sqrtRhoR = std::sqrt(rhoR);
  const tD vTilde = (sqrtRhoL * vL + sqrtRhoR * vR) / ( sqrtRhoL + sqrtRhoR);
  const tD hTilde = ((uL[2]+pL)/sqrtRhoL + (uR[2]+pR)/sqrtRhoR) / ( sqrtRhoL + sqrtRhoR);
  const tD aTilde = std::sqrt( (_cLaw.getGamma()-1.0) * (hTilde-0.5*uTilde*uTilde) );
  const tD sL = vTilde - aTilde;
  const tD sR = vTilde + aTilde;
  */

  /*
  // Toro's pressure based estimates: (10.59) -- (10.62)
  const tD gamma = _cLaw.getGamma();
  const tD rhoAve = 0.5*(rhoL+rhoR);
  const tD aAve   = 0.5*(aL  +aR  );
  const tD pPVRS  = 0.5*(pL+pR) - 0.5*(vR-vL)*rhoAve*aAve;
  tD qL = 1;
  if (pPVRS>pL)
    qL = std::sqrt(1+ (gamma+1.0)/(2*gamma)*(pPVRS/pL-1));
  const tD sL = vL - aL*qL;
  tD qR = 1;
  if (pPVRS>pR)
    qR = std::sqrt(1+ (gamma+1.0)/(2*gamma)*(pPVRS/pR-1));
  const tD sL = vR + aR*qR;
  */

  /* sStar: velocity of the contact wave (10.70) */
  const tD sStar = (pR - pL + rhoL*vL*(sL-vL) - rhoR*vR*(sR-vR)) /
                   (rhoL*(sL-vL) - rhoR*(sR-vR)) ;

  /* Compute fluxes */
  if (sL >=0){
    //F = Fleft
    F = _cLaw.getFprim(rhoL, vL, pL);
  }
  else if (sR <=0) {
    //F = Fright
    F = _cLaw.getFprim(rhoR, vR, pR);
  }
  // here sL<0 and sR>0 and use (10.74)
  else {
    t_u dStar;
    dStar[0]=0.;
    dStar[1]=1.;
    dStar[2]=sStar;
    if (sStar >=0){
      //F = FstarL
      F = sL*uL - _cLaw.getFprim(rhoL, vL, pL);
      F *= sStar;
      dStar *= (sL*(pL + rhoL*(sL-vL)*(sStar-vL)));
      F += dStar;
      F /= (sL-sStar);
    }
    else{
      //F = FstarR
      F = sR*uR - _cLaw.getFprim(rhoR, vR, pR);
      F *= sStar;
      dStar *= sR*(pR + rhoR*(sR-vR)*(sStar-vR));
      F += dStar;
      F /= (sR-sStar);
    }
  }
  lL = std::abs(sL);
  lR = std::abs(sR);
  return F;
}

hllcFlux<CLawLinear>::hllcFlux(CLawLinear & cLaw):
  _cLaw(cLaw)
  {}

/*! \brief The HLLC flux for the Linear transport equation
 *
 * The HLLC flux in this case is simply the upwind flux
 */
typename hllcFlux<CLawLinear>::t_u hllcFlux<CLawLinear>::getF(const t_u& uL/*! left state*/, const t_u& uR/*!right state*/, tC, tC, tD & sL/*!set to speed of left-going wave (or 0) in output*/, tD& sR/*!set to speed of right-going wave (or 0) in output*/)
{
  tD vel=_cLaw.getVel();
  if (vel>0){
    sL=0.;
    sR=vel;
    return vel*uL;
  }
  else{
    sL=-vel;
    sR=0.;
    return vel*uR;
  }
}

hllcFlux<CLawBurgers>::hllcFlux(CLawBurgers & cLaw):
  _cLaw(cLaw)
  {}

typename hllcFlux<CLawBurgers>::t_u hllcFlux<CLawBurgers>::getF(const t_u& uL/*! left state*/, const t_u& uR/*!right state*/, tC, tC, tD & lL/*!set to speed of left-going wave (or 0) in output*/, tD& lR/*!set to speed of right-going wave (or 0) in output*/)
{
  t_u F;
  tD sL = _cLaw.getVel(uL,0,0);
  tD sR = _cLaw.getVel(uR,0,0);

  /* Compute fluxes */
  if (sL >=0){
    //F = Fleft
    F = _cLaw.getF(uL,0,0);
  }
  else if (sR <=0) {
    //F = Fright
    F = _cLaw.getF(uR,0,0);
  }
  else {
    //FLL
    t_u FL = _cLaw.getF(uL,0,0);
    t_u FR = _cLaw.getF(uR,0,0);
    F = sR * FL - sL * FR + sL*sR * (uR - uL);
    F /= (sR - sL);
  }
  lL = std::abs(sL);
  lR = std::abs(sR);
  return F;
}

/** @} */
