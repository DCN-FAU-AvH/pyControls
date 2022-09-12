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

/** @addtogroup reconstructions
 *
 * @{
 */

/*! \file cweno.tpp
 *  \brief Implementation of cweno
 */

#ifndef CWENO_TPP
#define CWENO_TPP

#include <iomanip>

#include "cweno.hpp"
#include "udivdiff.hpp"
#include "polyTools.hpp"

/*! \brief Constructor from grid, accuracy and epsilon
 * 
 */
template <int M>
cweno<M>::cweno(Grid & grid, int accuracy, CWType type, tD epsilon, int alphaPower):
  RecBase<M>(grid),
  _udd(),
  _epsilon(epsilon),
  _type(type),
  _alphaPower(alphaPower)
  { setAccuracy(accuracy); }

/*! \brief Sets the accuracy of the reconstruction
 *
 * The accuracy has to be an odd integer. This in enforced here by
 * increasing the accuracy by 1 if an even integer is passed.
 * 
 * Resizes all internal storage vectors.
 */
template <int M>
void cweno<M>::setAccuracy(int accuracy){
  _accuracy = 2*(accuracy/2)+1; //makes sure it's odd...

  _nPolyLow = (_accuracy+1)/2; //Number of polys of low degree
  _polyLow.resize(_nPolyLow);

  //Resizes the storage of polyLow for polynomial coefficients
  _degLow = _nPolyLow-1;
  for (int k=0; k<_nPolyLow; ++k)
    _polyLow[k].resize(_degLow+1);
  
  _nPolyTot = _nPolyLow+1; //Number of total polys
  //Resizes the storage of optWeigths for optimal coefficients of cweno
  _optWeights.resize(_nPolyTot);
  
  setIdealWeights(_accuracy);
  
  //Resizes the storage of polyInds for indicators of the polys
  _polyInds.resize(_nPolyTot);
  
  //Resizes the storage of w nonlinear weights
  _w.resize(_nPolyTot);

  //Sets deafult value for alphaPower
  if (_alphaPower==0)
    setDefaultAlphaPower();

  _useIP0 = getDefaultIP0();
}

//! \brief Sets the value of the alphaPower parameter
template <int M>
void cweno<M>::setAlphaPower(int alphaPower)
  { _alphaPower = alphaPower; }

//! \brief Returns the value of the alphaPower parameter
template <int M>
int cweno<M>::getAlphaPower() const
  { return _alphaPower; }

//! \brief Sets the default value for the alphaPower parameter
template <int M>
void cweno<M>::setDefaultAlphaPower(){
  switch (_type){
  case CW:
    _alphaPower=2;
    break;
  case CWZdb:
    switch (_accuracy){
        case  3: _alphaPower=1; break;
        case  5: _alphaPower=2; break;
        default: _alphaPower=2;
    }
    break;
  case CWZ:
    switch (_accuracy){
      case  3: _alphaPower=1; break;
      case  5: _alphaPower=2; break;
      default: _alphaPower=2;
    }
    break;
  default:
    abort();
  }
}


//! \brief Returns the current accuracy of the reconstruction
template <int M>
int  cweno<M>::getAccuracy(){
  return _accuracy;
}

/*! \brief Sets the type of the type of CWENO reconstruction
 *
 */
template <int M>
void cweno<M>::setType(CWType type){
  _type=type;
}

//! \brief Returns the type of CWENO reconstruction
template <int M>
typename cweno<M>::CWType cweno<M>::getType() const{
  return _type;
}

/*! \brief Sets the epsilon of the reconstruction
 * 
 */
template <int M>
void cweno<M>::setEpsilon(tD epsilon){
  _epsilon=epsilon;
}

//! \brief Returns the value of epsilon
template <int M>
tD cweno<M>::getEpsilon() const{
  return _epsilon;
}

/*! \brief Sets the rule to compute I0
 *
 *  false (default): I0=I[Popt]
 *  true           : I0=I[P0]
 *
 * Note that true is the behaviour employed in early definitions of CWENO
 * and it still works better for order 3.
 */
template <int M>
void cweno<M>::setIP0(bool IP0){
  _useIP0 = IP0;
}

//! \brief Returns the rule to compute I0
template <int M>
bool cweno<M>::getIP0() const{
  return _useIP0;
}

//! \brief Returns the default value for _useIP0, depending on _type and _accuracy
template <int M>
bool cweno<M>::getDefaultIP0() const
{
  if (_type==CWZ)
    if (_accuracy<=3)
      return true;
    else
      return false;
  else
    return true;
}

/*! \brief Sets the optimal weights of the cweno reconstruction
 * 
 */
template <int M>
void cweno<M>::setIdealWeights(int accuracy){
  int k;
  
  _optWeights[0]=0.0;
  for(k=0;k<_nPolyLow/2;++k){
    _optWeights[1+k]=(1<<k);
    _optWeights[_nPolyLow-k]=_optWeights[1+k];
    _optWeights[0] += 2*_optWeights[1+k];
  }
  if (_nPolyLow%2){
    _optWeights[1+k]=(1<<k);
    _optWeights[0] += _optWeights[1+k];
  }
  tD d_tot=4*_optWeights[0];
  _optWeights[0]*=3.;
  for(k=0;k<_nPolyTot;++k)
    _optWeights[k] /= d_tot;
}

/*! \brief Computes all the reconstruction polynomials
 * 
 */
template <int M>
void cweno<M>::compute(const DofVec<M> * u, int cellStart, int cellEnd){
  setUAvgAndResize(u);
  _udd.compute(*u,_accuracy-1);

  int r = _accuracy/2; //shift for central polynomial

  if ( (_polyHigh.size() < u->size() ) ||
       (_polyHigh[0].size() < _accuracy)                       )
    setAccuracy(_accuracy);

  //Loop on the grid
  for (int k=cellStart; k<cellEnd; ++k){
    //Compute the optimal polynomial on the central stencil
    _udd.InterpPoly(k, _accuracy-1, r, _polyHigh[k].data() );
    if (! _useIP0)
      _polyInds[0]=computeIndicator(_accuracy-1, _polyHigh[k].data());
    //Compute low order polynomials and indicators
    for(int i=0;i<_nPolyLow;++i){
      _udd.InterpPoly(k, _degLow, i, _polyLow[i].data() );
      _polyInds[i+1]=computeIndicator(_degLow, _polyLow[i].data());
    }
    //Compute P0 and the indicator
    for(int j=0;j<_accuracy;++j)
      _polyHigh[k][j] = _polyHigh[k][j] / _optWeights[0];
    for(int i=0;i<_nPolyLow;++i)
      for(int j=0;j<=_degLow;++j)
        _polyHigh[k][j]=_polyHigh[k][j]-(_optWeights[i+1]/_optWeights[0])*_polyLow[i][j];
    if (_useIP0)
      _polyInds[0]=computeIndicator(_accuracy-1, _polyHigh[k].data());
    //Compute nonlinear weights
    switch (_type){
      case CW:
        for(int i=0;i<_nPolyTot;++i)//{
          _w[i]= _optWeights[i] / eWisePow(_polyInds[i]+_epsilon, _alphaPower);
      break;
      case CWZdb:
        switch (_accuracy){
          case 3:
            _tau = eWiseAbs(_polyInds[1] - _polyInds[2]);
            /* simple choice for tau does not exist */
          break;
          case 5:
            _tau = eWiseAbs(_polyInds[3] - _polyInds[1]);
            /* simple choice for tau is the same as the above optimal one */
          break;
          case 7:
            _tau = eWiseAbs(_polyInds[1] + tD(3)*_polyInds[2] - tD(3)*_polyInds[3] - _polyInds[4]);
            /* tau = eWiseAbs(_polyInds[1] - _polyInds[2] - _polyInds[3] + _polyInds[4]); /* simple choice for tau */
          break;
          case 9:
            _tau = eWiseAbs(_polyInds[1] + tD(2)*_polyInds[2] - tD(6)*_polyInds[3] + tD(2)*_polyInds[4] + _polyInds[5]);
            /* tau = eWiseAbs(_polyInds[_nPolyTot-1] - _polyInds[1]); /* simple choice for tau */
          break;
          case 11:
            _tau = eWiseAbs(_polyInds[1] + _polyInds[2] - tD(8)*_polyInds[3] + tD(8)*_polyInds[4] - _polyInds[5] - _polyInds[6]);
            /* tau = eWiseAbs(_polyInds[1] - _polyInds[2] - _polyInds[3] + _polyInds[4]); /* simple choice for tau */
          break;
          case 13:
            _tau = eWiseAbs(_polyInds[1] + tD(36)*_polyInds[2] + tD(135)*_polyInds[3] - tD(135)*_polyInds[5] - tD(36)*_polyInds[6] - _polyInds[7]);
            /* tau = eWiseAbs(_polyInds[_nPolyTot-1] - _polyInds[1]); /* simple choice for tau */
          break;
          default:
            std::cout << "Value of tau not known in Don-Borges(2013) paper for CWENOZ" << _accuracy << ". Aborting." << std::endl;
            abort();
        }
        for(int i=0;i<_nPolyTot;++i)
          _w[i] = _optWeights[i] * ( eWisePow( eWiseDivide( _tau, _polyInds[i]+_epsilon) , _alphaPower) + tD(1.));
#ifdef DEBUG_REC
        std::cerr << _tau[0] << ", ";
#endif
      break;
      case CWZ:
        switch (_accuracy){
          case 3:
            _tau = eWiseAbs(_polyInds[1] + _polyInds[2] - tD(2)*_polyInds[0]);
          break;
          case 5:
            _tau = eWiseAbs(_polyInds[1] + tD(4)*_polyInds[2] + _polyInds[3] - tD(6)*_polyInds[0]);
          break;
          default:
            std::cout << "Value of tau not computed for CWENOZ yet" << _accuracy << ". Aborting." << std::endl;
            abort();
        }
        for(int i=0;i<_nPolyTot;++i)
          _w[i] = _optWeights[i] * ( eWisePow( eWiseDivide( _tau, _polyInds[i]+_epsilon) , _alphaPower) + tD(1.));
#ifdef DEBUG_REC
        std::cerr << _tau[0] << ", ";
#endif
      break;
      default:
        std::cout << "Unknown CWENO type " << _type << ". Aborting." << std::endl;
        abort();
    }
    t_u sum_w(0.0);
    for(int i=0;i<_nPolyTot;++i)
      sum_w+=_w[i];
    for(int i=0;i<_nPolyTot;++i)
      _w[i]/=sum_w;
#ifdef DEBUG_REC
    for(int i=0;i<_nPolyTot;++i)
      std::cerr << std::abs(_w[i][0]-_optWeights[i]) << ", ";
#endif
    //Compute the reconstruction polynomial
    for(int j=0;j<_accuracy;++j)
      _polyHigh[k][j]=eWiseProduct(_w[0],_polyHigh[k][j]);
    for(int i=0;i<_nPolyLow;++i)
      for(int j=0;j<=_degLow;++j)
        _polyHigh[k][j]=_polyHigh[k][j]+eWiseProduct(_w[i+1],_polyLow[i][j]);
  }
}

/*! \brief evaluates the reconstruction inside the cell.
 * 
 * pos should be between -0.5 and 0.5.
*/
template <int M>
typename cweno<M>::t_u cweno<M>::eval(int k, tC pos){
  return evalPoly(_accuracy-1, pos, _polyHigh[k].data() );
}

template <int M>
const typename cweno<M>::t_u * cweno<M>::getCoeffRec(int k) const{
  return _polyHigh[k].data();
}

template <int M>
void cweno<M>::overrideCoeffRec(int k, const t_u* coeff){
  for (int d=0; d<_accuracy+1; ++d)
    _polyHigh[k][d] = coeff[d];
}

//! \brief Returns the number of ghost cells used by this reconstruction
template <int M>
int cweno<M>::needsGhosts() const{
  return _accuracy/2;
}

//! \brief Returns the order of accuracy
template <int M>
int cweno<M>::getOrder() const{
  return _accuracy;
}

//! \brief Resize storage for reconstruction polynomial
template <int M>
void cweno<M>::setUAvgAndResize(const DofVec<M> * u){
  if ( (_polyHigh.size() < u->size() ) || (_polyHigh[0].size() < _accuracy) ) {
    _polyHigh.resize(u->size());
    for (int k=0; k<u->size(); ++k)
      _polyHigh[k].resize(_accuracy+1);
  }
}

//! Prints info on the reconstruction
template <int M>
void cweno<M>::print() const {
  switch (_type){
  case CW:
    std::cout << "CWENO "; break;
  case CWZdb:
    std::cout << "CWENOZdb "; break;
  case CWZ:
    std::cout << "CWENOZ "; break;
  default: abort();
  }
  std::cout << " reconstruction of order " << _accuracy;
  std::cout << ", using alphaPower=" << _alphaPower
            << ", epsilon=" << std::scientific << std::setprecision(3) << _epsilon;
  if ((_type==CWZdb) || (_type==CWZ))
    std::cout << " and I0=I[" << (_useIP0?"P0":"Popt") << "]";
  std::cout << std::endl;
};


#endif

/** @} */

