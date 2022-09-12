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

/** @addtogroup utils
 *
 * @{
 */

/*! \file parseOptions.cpp
 *  \brief Implementation of Options class
 */

#include "parseOptions.hpp"
#include <iostream>

/*! \brief Constructor: sets the known options */
Options::Options():
  desc("Allowed options")
  {
    desc.add_options()
      ("help", "produce help message")
      ("N", po::value<int>(&_Ncells)->default_value(25), "N of cells")
      ("cfl", po::value<double>(&_cfl)->default_value(0.45), "N of cells")
      ("pb", po::value<std::string>(&_pbString)->default_value("default"), "Problem name")
      ("quad", po::value<int>(&_quadOrder)->default_value(0), "quadrature order for the initial data")
      ("rec", po::value<std::string>(&_recName)->default_value("const"), "Reconstruction")
      ("recalphapow", po::value<int>(&_recAlphaPower)->default_value(0), "power parameter in reconstructions")
      ("recepspow", po::value<int>(&_recEpsPower)->default_value(0), "power for choosing epsilon in reconstructions (power of dx if positive, power of 10 if negative)")
      ("recforceIP0", po::bool_switch(&_forceIP0)->default_value(false), "Use I0=I[P0]")
      ("recforceIPopt", po::bool_switch(&_forceIPopt)->default_value(false), "Use I0=I[Popt]")
      ("charProj", po::bool_switch(&_charProj)->default_value(false), "Use characteristic projection")
//#ifdef SWE
      ("WBhydro", po::bool_switch(&_WBhydro)->default_value(false), "Use hydrostatic well-balancing")
//#endif
      ("flux", po::value<std::string>(&_fluxName)->default_value("llf"), "Numerical Flux")
      ("stepper", po::value<std::string>(&_timeStepperName)->default_value("euler"), "Time Stepper")
      ("tfin", po::value<double>(&_tFinal), "override default final integration time");
  }

/*! \brief Parses the command line, stores the values of the known
 * options and prints help if required
 */
void Options::parseCmdLine(int argc, const char* const argv[]){
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (_forceIP0 && _forceIPopt){
    std::cout << "Only one of --recforceIP0 and --recforceIPopt may be given!" << std::endl;
    abort();
  }
  _forcedI0 = _forceIP0 | _forceIPopt;
  _useIP0 = _forceIP0 | (!_forceIPopt);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(0);
  }
}

//! \brief Returns the value of the --cfl option
tD Options::getCFL()
  { return _cfl;}

//! \brief Returns the value of the --N option
int Options::getNcells()
  { return _Ncells;}

//! \brief Returns the value of the --quad option
int Options::getQuadOrder()
  { return _quadOrder;}

//! \brief Returns the value of the --tfin option
bool Options::getFinalTime(tC &tFin){
  if (vm.count("tfin")){
    tFin = _tFinal;
    return true;
  }
  else
    return false;
}

//! \brief Returns the value of the --pb option
std::string Options::getPBString()
{ return _pbString;}

//! \brief Returns the value of the --rec option
std::string Options::getRecName()
{ return _recName;}

//! \brief Returns the value of the --recalphapow option
int Options::getRecAlphaPower()
  { return _recAlphaPower;}

//! \brief Returns the value of the --recepspow option
int Options::getRecEpsPower()
  { return _recEpsPower;}

//! \brief Returns true if the command option -recuseIP0 was given
bool Options::forcedI0()
{
  return _forcedI0;
}

//! \brief Returns the value of the option -recuseIP0. Warning: check givenuseIP0() first!
bool Options::getuseIP0()
{
#ifdef DEBUG
  if (!getuseIP0())
    abort();
#endif
  return _useIP0;
}

//! \brief Returns the value of the --charProj option
bool Options::getCharProj(){
  return _charProj;
}

//#ifdef SWE
//! \brief Returns the value of the --WBhydro option
bool Options::getWBhydro(){
  return _WBhydro;
}
//#endif

//! \brief Returns the value of the --flux option
std::string Options::getFluxName()
{ return _fluxName;}

//! \brief Returns the value of the --stepper option
std::string Options::getTimeStepperName()
{ return _timeStepperName;}
