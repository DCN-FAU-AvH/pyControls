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

/** @addtogroup grid
 *
 * @{
 */

/*! \file grid.cpp
 *  \brief Definition of Grid
 */

/** @} */

#include "grid.hpp"

/*! \brief Constructor
 *
 * Input: x0,x1 left and right boundary
 *        N  number of physical cells
 *        gh number of ghost cells (optional: defaults to 0)
 */
Grid::Grid(tC x0, tC x1, int N, int gh):
    _N(N), _gh(gh), _xLeft(x0), _xRight(x1), _dx((x1-x0)/N)
    {}

//! \brief Sets the number of ghosts
int Grid::setGhosts(int gh){
  _gh=gh;
}

//! \brief Number of ghosts
int Grid::getGhosts() const
  {return _gh;}
//! \brief Number of physical cells
int Grid::getPhysicalSize() const
  {return _N;}
//! \brief Number of cells + one ghost per side
int Grid::getOneGhostSize() const
  {return (_N+2);}
//! \brief Number of cells (including ghosts)
int Grid::getFullSize() const
  {return (_N+2*_gh);}
//! \brief Number of faces (both internal and domain boundary)
int Grid::getFluxesSize() const
  {return (_N+1);}

bool Grid::isGhost(int cellId) const
  {
    if ((cellId<_gh) || (cellId>_gh+_N) )
      return true;
    return false;
  }
bool Grid::isPhysical(int cellId) const
  {
    if ((cellId<_gh) || (cellId>_gh+_N) )
      return false;
    return true;
  }

//! \brief id first ghost on the left
int Grid::beginOneGhost() const
  {return _gh-1;}
//! \brief id one past the first ghost on the right
int Grid::endOneGhost() const
  {return _gh+_N+1;}

//! \brief id of first physical cell on the left
int Grid::beginPhysical() const
  {return _gh;}
//! \brief id of one past the last physical cell on the right
int Grid::endPhysical() const
  {return _gh+_N;}

//! \brief id of first cell on the left
int Grid::beginFull() const
  {return 0;}
//! \brief id of one past the last cell on the right
int Grid::endFull() const
  {return _N+2*_gh;}

//! \brief id of first ghost cell on the right
int Grid::beginRightGhosts() const
  {return _gh+_N;}
//! \brief id of one past the last ghost cell on the right
int Grid::endRightGhosts() const
  {return _N+2*_gh;}

/*! \brief id of first ghost cell on the left
 *
 * Use with counter-- increments!
 */
int Grid::rbeginLeftGhosts() const
  {return _gh-1;} //first ghost cell on the left
/*! \brief id of one before the last ghost cell on the left
 *
 * Use with counter-- increments!
 */
int Grid::rendLeftGhosts() const
  {return -1;} //one before the last ghost cell on the left

//! \brief id of first face
int Grid::beginFaces() const
  {return 0;} //first face
//! \brief id of one past the last face
int Grid::endFaces() const
  {return _N+1;} //one past the last face

//! \brief center of cell
tC Grid::getXCenter(int cellId) const
  {return (_xLeft-_gh*_dx+(cellId+0.5)*_dx);}
//! \brief size of cells
tC Grid::getDx() const
  {return _dx;}

//! \brief position of face
tC Grid::getXFace(int faceId) const
  {return (_xLeft+faceId*_dx);}
//! \brief id of cell on the left of a given interface
int Grid::leftOfFace(int faceId) const
  {return beginOneGhost()+faceId;}
//! \brief id of cell on the right of a given interface
int Grid::rightOfFace(int faceId) const
  {return beginPhysical()+faceId;}
//! \brief id of face on the left of a given cell
int Grid::leftFace(int cellId) const
  {return beginOneGhost()+cellId;}
//! \brief id of face on the left of a given cell
int Grid::rightFace(int cellId) const
  {return beginPhysical()+cellId;}

//! \brief set boundary condition for left boundary
void Grid::setLeftBC(BCType lBC)
  {leftBC = lBC;}
//! \brief get boundary condition for left boundary
Grid::BCType Grid::getLeftBC() const
  {return leftBC;}
//! \brief set boundary condition for right boundary
void Grid::setRightBC(BCType rBC)
  {rightBC = rBC;}
//! \brief get boundary condition for right boundary
Grid::BCType Grid::getRightBC() const
  {return rightBC;}
