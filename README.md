# pyControls

The pyControls package provides a framework for simulations of gas networks with arbirary geometry and edge dynamics.
This package provides a high-order solver for the isoeuler equations and a reader for the GasLib scenarios. Furthermore,
the network is easily plotted either in 2D or 3D manner.

<p align="center">
<img src="https://github.com/DCN-FAU-AvH/pyControls/blob/main/Icon.png" width="70%" height="70%" >
</p>

## Quick installation

In order to install pyControls, run

```
pip install pyControls
```
(Todo: the package will be defined when the 0.1 version of the code is finished)

Install python prerequisites:
```
pip3 install networkx pint numpy scipy matplotlib cppyy
```

~~At the moment:~~
```
cd lib/claw1darena
mkdir build
cd build
cmake ..
make -j
```
~~Then run the `src/network.py` file for the diamond network example. (Set the correct pathes in cw.py and network.py)~~

Use the crossout version if you're at least a little familiar with C++. Otherwise use the local Lax-Friedrich Flux version. After installing the python prerequisistes, run `src/main_diamond.py`.



## Usage 

### Simple network of two nodes and one edge

Defining a simple network consisting of two nodes is done as following. First import the pyControls package together with numpy and math, then define one edge. The edge dynamics requires setting the parameters of the underlying physical system (in this case the isoeuler equations), geometrical information and the cfl number.

```python
import math
import numpy as np
import matplotlib.pyplot as plt
import network as nw

# import iseuler stuff
from isoeuler import llf, pde, nodes

pipe_id = 0
L = 1000 # length of the pipe in m
dx = 1 # cell width
cfl = 0.9
alpha = 1.0 #
u0 = lambda x : x > 500 ? np.array([1, 0]) : np.array([0.5, 0]) # initial condition Riemann problem
p1 = nw.Pipe(pipe_id, llf.iso_llf(pipe_id,
                                          math.floor(L/dx),
                                          L, cfl,
                                          alpha, gamma,
                                          u0))
```
The `nw.Pipe` object is constructed by passing an id and a dynamics object, for the definition of latter one see the development wiki. In this case we use a first order FV scheme with local Lax-Friedrich flux for the isoeuler equations.
The isoeuler system requires the alpha and beta parameters for the equations of state p(rho) = alpha * rho^gamma.
Each edge is defined on the interval [0, L]. 

Next we need to define two nodes. In this case we've used Boundary nodes for the isoeuler system.

```python
n1 = nodes.Boundary(1, # id
              lambda u, t : u,  # boundary law
              {pipe_id : np.array([-1,0])}, # vector nu
              {pipe_id : 0}) # direction of 
n2 = nodes.Boundary(2, 
              lambda u, t : u, 
              {pipe_id : np.array([1,0])}, 
              {pipe_id : 1})
```
The boundary law may be set arbitrary by passing a function depending on the state vector of the adjacent
pipe `u` and time `t`.  The vector `nu` gives the direction of the pipe relative to the node and its 2-norm
encodes the thickness of the adjacent pipe. Finally, the Boundary object requires an indication which end
of the pipe is connected to the node, if `pipe_id:0` the pipe with id `pipe_id` is connected with its left-hand side boundary to the
node and setting `pipe_id:1` indicates that the pipe is connected with the right-hand side boundary.

Finally, we collect the edges and nodes in dictionaries and construct a network object.
```python
nw1 = nw.Network([n1, n2], {(0,1) : p1})
```
The edges of a network are stored in dictionary where a pair of node ids is mapped to the pipe object.

To run a calculation we define a `while` loop for the time and call the `nw1.dvance` repeatedly.

```python
while t < Tfinal:
    nw1.advance()
```

For the plots of the network the matplotlib library is employed. Defining a figure `fig` and axis `axs` objects, the network can be plotted as a 3D-lines plot using `nw.plot_3d(nw1, fig, axs)` or 2D-network using `nw.plot(nw1, fig, axs)`.

Adding more nodes and edges to the above defined `nw` object allows calculating arbitrary gas networks.

### GasLib reader

The pyControls package provides a reader for the GasLib XML format. The examples from the GasLib library can
directly be imported as zip archives and creates a network object using the isoeuler dynamics. 
```python
nw = pyControls.Gaslib.gaslib_to_isoeuler('path_to_gaslib_archive.zip')
```
The reader can be extended to other dynamics, see the development documentation.

## Development documentation

In order to add a new PDE system together with its node and pipe dynamics it is recommended to follow the same pattern as for the isoeuler equations.
Make a folder for your new system and therin place three modules for the pde, pipe dynamics and nodes, respectively.

### Adding a node type
Derive a class from the `Coupling_node` base class and implement the `calc` function. The `calc` function receives all connected pipe objects, calculates the coupled states and sets the boundary conditions in each pipe solver.

### Adding pipe dynamics
Implementing a first order FV scheme is easily done by deriving from the `fv` class and implementing the numerical flux `nf` and the source term discretization `source` methods. Following the structure of the `fv` class one can implement an arbitrary solver that should be able to advance, calculate timestep sizes, provide an interface to the solution and set the boundary conditions.

<!---### Extending the GasLib reader-->

<!---## todos-->
