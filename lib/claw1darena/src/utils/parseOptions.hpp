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

/** @defgroup utils Sundry useful tools
 *
 * @{
 */

/*! \file parseOptions.hpp
 *  \brief Declaration of Options class
 */

#ifndef PARSEOPTIONS_HH
#define PARSEOPTIONS_HH

#include "../config.h"

#include <boost/program_options.hpp>
#include <string>

namespace po = boost::program_options;

/*! \brief Parsing of command line options
 *
 * This class has a method to parse command line options and stores
 * internally the values found.
 *
 * In order to add a new command line option, it should be added to
 * thre desc variable in the constructor.
 */
class Options{
public:
  Options();
  void parseCmdLine(int argc, const char* const argv[]);
  tD getCFL();
  int getNcells();
  int getQuadOrder();
  bool getFinalTime(tC &tFin);
  std::string getPBString();
  std::string getRecName();
  int getRecAlphaPower();
  int getRecEpsPower();
  bool getCharProj();
  bool getuseIP0();
  bool forcedI0();
//#ifdef SWE
  bool getWBhydro();
//#endif
  std::string getFluxName();
  std::string getTimeStepperName();

private:
  double _cfl;    //!< storage for cfl number (as double, due to Boost requirement)
  double _tFinal; //!< storage for final time (as double, due to Boost requirement)
  int _Ncells;    //!< storage for the number of cells
  int _quadOrder; //!< storage for the quadrature order
  std::string _pbString; //!< storage for the problem name
  std::string _recName;  //!< storage for the reconstruction name
  int _recAlphaPower; //!< storage for the power parameter in reconstructions
  int _recEpsPower; //!< storage for the power parameter in reconstructions (epsilon=dx^_recEpsPower if positive, epsilon=10^(_recEpsPower) if negative)
  bool _forceIP0;   //!< force use of I0=I[P0]
  bool _forceIPopt; //!< force use of I0=I[Popt]
  bool _forcedI0; //!< was I0 specified in the options?
  bool _useIP0;   //!< use I0=I[P0]? (instead of I[Popt])
  bool _charProj; //!< use characteristic projection?
//#ifdef SWE
  bool _WBhydro;    //!< use well-balanced hydrostatic method
//#endif
  std::string _fluxName; //!< storage for the numerical flux name
  std::string _timeStepperName; //!< storage for the timestepper name

  po::options_description desc; //!< known options
  po::variables_map vm;
};

#endif
