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

/** @defgroup grid Representation of a grid and related tools
 *
 * @{
 */

/*! \file grid.hpp
 *  \brief Declaration of Grid
 */

#ifndef GRID_HH
#define GRID_HH

#include "../config.h"

/*! \brief Lightweight class to represent a finite volume uniform grid in 1d.
 *
 * A grid on an interval [x0,x1] is composed of N cells and has gh
 * ghost cells per side.
 *
 * The user is assumed to store cell averages data in a DoFVec with
 * indices starting from 0.
 *
 * Data attached to faces (e.g. fluxes) should also be stored in a DoFVec
 * with indices starting from 0.
 *
 * Methods called ***Size return the number of faces, cells, physical
 * cells, etc. These numbers should be used to allocate the
 * appropriate amount of storage.
 *
 * Several pairs of Begin/End methods are provided to ease the task of
 * looping over several meaningful sets of cells. Also, convenience
 * functions that return the cell indices on right/left of a given face
 * and the indices of the faces on the right/left of a given cell are
 * provided.
 *
 * IMPORTANT NOTE: the number of ghosts can be changed after grid
 * creation by the method setGhosts() and of course the number of ghosts
 * influence the return values of all the begin/end*** and ***Size
 * methods. If you call setGhosts after using the output of any of the
 * size or begin/end methods, you're looking for troubles...
 *
 */
class Grid{
public:
  Grid(tC x0, tC x1, int N, int gh=0);

  /*! \brief Enumeration of boundary conditions
   *
   * Note: a CLaW class will recognize and implement only those that
   * are meaningful in that context.
   */
  enum BCType {
    BC_PERIODIC  = 0, /*!< Periodic boundary */
    BC_FREEFLOW  = 1, /*!< Free-flow boundary */
    BC_DIRICHLET = 2, /*!< Dichlet boundary */
    BC_WALL      = 3  /*!< Wall boundary */
  };

  int setGhosts(int gh);
  int getGhosts() const;
  int getPhysicalSize() const;
  int getOneGhostSize() const;
  int getFullSize() const;
  int getFluxesSize() const;

  bool isGhost(int cellId) const;
  bool isPhysical(int cellId) const;

  int beginOneGhost() const;
  int endOneGhost() const;

  int beginPhysical() const;
  int endPhysical() const;

  int beginFull() const;
  int endFull() const;

  int beginRightGhosts() const;
  int endRightGhosts() const;

  int rbeginLeftGhosts() const;
  int rendLeftGhosts() const;

  int beginFaces() const;
  int endFaces() const;

  tC getXCenter(int cellId) const;
  tC getDx() const;

  tC getXFace(int faceId) const;
  int leftOfFace(int faceId) const;
  int rightOfFace(int faceId) const;
  int leftFace(int cellId) const;
  int rightFace(int cellId) const;

  void setLeftBC(BCType lBC);
  BCType getLeftBC() const;
  void setRightBC(BCType rBC);
  BCType getRightBC() const;

private:
  const int _N;     //!< number of internal cells
  int _gh;          //!< number of ghost cells per side
  const tC _xLeft;  //!< location of left domain boundary
  const tC _xRight; //!< location of right domain boundary
  const tC _dx;     //!< cell sizes
  BCType leftBC;    //!< left  boundary type
  BCType rightBC;   //!< right boundary type
};

#endif

/** @} */
