#!/usr/bin/python

#  This file is part of claw1dArena, a software library to discretize
#  one-dimensional conservation laws.
#
#  Copyright (C) 2017 Matteo Semplice
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#  For details on the licence, see the COPYING file at the root of the tree.

import matplotlib.pyplot as plt
import numpy as np

x0, h0, q0, z0 = np.loadtxt('initial.csv', delimiter=' ', unpack=True)
x1, h1, q1, z1 = np.loadtxt('final.csv', delimiter=' ', unpack=True)

plt.subplot(2,1,1)
plt.plot(x0,h0+z0, label='initial')
plt.plot(x1,h1+z1, label='final')
plt.plot(x0,z0, label='bottom')
plt.xlabel('x')
plt.title('height')
plt.grid()

v0=q0/h0;
v1=q1/h1;
plt.subplot(2,1,2)
plt.plot(x0,v0, label='initial')
plt.plot(x1,v1, label='final')
plt.xlabel('x')
plt.title('Velocity')
plt.grid()

plt.legend()
plt.show()
