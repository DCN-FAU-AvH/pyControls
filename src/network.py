from cw import * 
import networkx as nx
from isoeuler import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize


class Pipe():
    """ Pipe container  
    
    Attributes
    ----------
    id : int/str
     Unique id of the pipe
    solver : class
     solver for a single pipe dynamics, see solver.py for the specifications of
     the required attributes and methods of the solver classes and num.py
     for an example FV implementation
    """
    def __init__(self, id, solver):
        """ 
        Parameters
        ----------
        id : int/str
         Unique id of the pipe
        solver : class
         solver for a single pipe dynamics
        """
        self.id = id
        self.solver = solver

            
class Network():
    """ Network of pipes, nodes and boundary leaves

    Basically a wrapper around a networkx network

    Attributes
    ----------
    t : float
     current physical time in the network
    G : networkx.Graph
     networkx object consisting of pipes, nodes and boundary leaves

    Methods
    -------
    min_dt()
     calculates the minimal timestep over all pipes in the network, i.e.
     provides the correct timestep according to the cfl condition across
     the network
    
    advance()
     calculates the node and boundary conditions and sets the boundary
     conditions in the pipes. Then advances each pipe by the timestep
     given by min_dt
    """

    def __init__(self, nodes, edges_map):
        """

        Parameters
        ----------
        nodes : list of Coupling_node type objects
        edges_map :  dict (pair of ids of nodes) -> solver
         the order of the nodes in the map determines the direction of the cells
         in the pipe relative to the node
        """
        self.t = 0.0
        # todo: check that nodes is a list and edges is a dict
        self.G = nx.Graph()
        self.G.add_nodes_from(nodes)
        # add nodes
        for i in range(0, len(nodes)):
            self.G.add_node(nodes[i], id=i)
        # add edges
        for ids_pair, p in edges_map.items():
            self.G.add_edge(nodes[ids_pair[0]], nodes[ids_pair[1]], pipe=p)

    # methods
            
    def min_dt(self):
        """
        calculates the minimal timestep over all pipes in the network, i.e.
        provides the correct timestep according to the cfl condition across
        the network
        """
        dt = float('inf')
        for n1, n2, data in self.G.edges().data():
            dt = min(dt, data['pipe'].solver.get_dt())
        return dt

    def advance(self):
        """
        calculates the node and boundary conditions and sets the boundary
        conditions in the pipes. Then advances each pipe by the timestep
        given by min_dt
        """
        # todo: add dtmax
        # calculate all nodes conditions
        for n in self.G.nodes():
            adj_pipes=[self.G.edges[e] for e in self.G.edges(n)]
            n.calc(adj_pipes, self.t)
        # advance each edge
        dt = self.min_dt()
        for n1, n2, data in self.G.edges().data():
            data['pipe'].solver.set_dt(dt)
            data['pipe'].solver.advance(self.t)

        self.t += dt 

# todo: divide in smaller fcts,
# todo: use animation object of matplotlib! (check if FuncAnimation or TimedAnimation is better!)
# todo: add a function mapping the conservative variables to the quantity to be plotted
# todo: add to netork class and call the attributes directly
def plot_nw(nw, fig, axs):
    """ Plots the network as 2D colorbars 
    Parameters
    ----------
    nw : Network
     Network object, see above
    fig : Figure
     matplotlibs figure object to plot the nw to
    axs : Axis
     matplotlibs axis object to plot the nw to
    """
    edges_visited = set()
    pos = {}

    # move to sep. fct
    min_val =  0.3#float('inf') 
    max_val =  1.0#float('-inf')
    # todo: could be more efficient, maybe write some kind of an iterator for all edges in network?
    for n, nbrs in nw.G.adj.items():
        for nbr, eattr in nbrs.items():
                solver = eattr['pipe'].solver
                vals = np.array([solver.U[i][0] for i in  solver.physical_range])
#                min_val = min(min_val, min(vals))
#                max_val = max(max_val, max(vals))

                
    all_segments=np.empty((0,2,2))    
    all_vals = np.array([])
    # loop nodes and collect segments and vals
    for n, nbrs in nw.G.adj.items():
        for nbr, eattr in nbrs.items():
            if not((n, nbr) in edges_visited) and not((nbr, n) in edges_visited):
                edges_visited.add((n, nbr))
                solver = eattr['pipe'].solver
                nu = n.nus[eattr['pipe'].id]
                if n not in pos:
                    pos[n] = np.array([0.0, 0.0])
                if nbr not in pos:
                    pos[nbr] = pos[n] + nu/(np.linalg.norm(nu))*solver.L
                    
                # todo: put the following in a separate fct
                # todo replace by np.array slice
                vals = np.array([solver.U[0][i] for i in  solver.physical_range])
                x = np.linspace(pos[n][0], pos[nbr][0], len(vals) )
                y = np.linspace(pos[n][1], pos[nbr][1], len(vals) )
                # the following is from here https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/multicolored_line.html
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                all_segments = np.append(all_segments, segments, axis = 0)
                if n.dirs[solver.id] == 0:
                    all_vals = np.append(all_vals, vals)
                else:
                    all_vals = np.append(all_vals, np.flip(vals))

    norm = plt.Normalize(min_val, max_val)                
    lc = LineCollection(all_segments, cmap='coolwarm', norm=norm)
    lc.set_array(all_vals)
    # Set the values used for colormapping in the correct direction
    lc.set_linewidth(7.5*np.linalg.norm(nu))
    line = axs.add_collection(lc)

    axs.set_xlim(-2, 50)
    axs.set_ylim(-20, 20)
    return line, min_val, max_val

# todo: divide in smaller fcts,
# todo: use animation object of matplotlib! (check if FuncAnimation or TimedAnimation is better!)
# todo: add a function mapping the conservative variables to the quantity to be plotted
# todo: add to network class and call the attributes directly
# todo: unify with the plot_nw function
def plot_nw3d(nw, fig, axs):
    """ Plots the network as 3D colored lines 
    Parameters
    ----------
    nw : Network
     Network object, see above
    fig : Figure
     matplotlibs figure object to plot the nw to
    axs : Axis
     matplotlibs axis object to plot the nw to
    """
    edges_visited = set()
    pos = {}

    # move to sep. fct
    min_val =  0.3#float('inf') 
    max_val =  0.8#float('-inf')
    # todo: could be more efficient, maybe write some kind of an iterator for all edges in network?
    for n, nbrs in nw.G.adj.items():
        for nbr, eattr in nbrs.items():
            solver = eattr['pipe'].solver
            vals = np.array([solver.U[0, i] for i in  solver.physical_range])
            #                min_val = min(min_val, min(vals))
            #                max_val = max(max_val, max(vals))
                
    # loop nodes and collect segments and vals
    for n, nbrs in nw.G.adj.items():
        for nbr, eattr in nbrs.items():
            if not((n, nbr) in edges_visited) and not((nbr, n) in edges_visited):
                edges_visited.add((n, nbr))
                solver = eattr['pipe'].solver
                nu = n.nus[eattr['pipe'].id]
                if n not in pos:
                    pos[n] = np.array([0.0, 0.0])
                if nbr not in pos:
                    pos[nbr] = pos[n] + nu/(np.linalg.norm(nu))*solver.L
                    
                # todo: put the following in a separate fct
                # todo replace by np.array slice
                print(solver.U.shape, solver.physical_range)
                vals = np.array([solver.U[0][i] for i in  solver.physical_range])
                x = np.linspace(pos[n][0], pos[nbr][0], len(vals) )
                y = np.linspace(pos[n][1], pos[nbr][1], len(vals) )
                # the following is from here https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/multicolored_line.html
                if n.dirs[solver.id] == 0:
                    axs.plot(x[0:-2],y[0:-2],vals[0:-2])
                else:
                    axs.plot(x[0:-2],y[0:-2],np.flip(vals[0:-2]))
                # add geoemtry on the 'floor'
                axs.plot(x[0:-2],y[0:-2], np.zeros(len(x[0:-2])), 'k-')
                    
    axs.set_xlim(-2, 50)
    axs.set_ylim(-20, 20)
    axs.set_zlim(0, max_val)

    axs.set_xlabel('x')
    axs.set_ylabel('y')
    axs.set_zlabel('Density')


# todo:  star shaped network test, move to tests
# u0 = lambda x : [1.0, 0.0, 0.0] if x > 5.0 else [0.125, 0.0, 0.0]
# u1 = lambda x : [0.125, 0.0, 0.0]

# test_edges = {
#     (0,1) : Pipe(0, cweno3(0, nr_cells, L, cfl, alpha, gamma, u0)),
#     (0,2) : Pipe(1, cweno3(1, nr_cells, L, cfl, alpha, gamma, u1)),
#     (0,3) : Pipe(2, cweno3(2, nr_cells, L, cfl, alpha, gamma, u1)),
#     (4,0) : Pipe(3, cweno3(3, nr_cells, L, cfl, alpha, gamma, u1))
# }



# n1 = Junction("0", sys, {0 : np.array([1.0, 0.5]), 1: np.array([0.0, 1.0]), 2: np.array([0.0, -1.0]), 3: np.array([-1.0, 0.0])}, {0:0,1:0,2:0,3:0})
# #n1 = Junction({0 : np.array([1.0, 0.0]), 1: np.array([0.0, 1.0])}, {0:0,1:0,2:0,3:0})
# n2 = Boundary("1", sys, lambda u,t : u, {0: np.array([-1.0, 0.5])}, {0:1})
# n3 = Boundary("2", sys, lambda u,t : u, {1: np.array([0.0, -1.0])}, {1:1})
# n4 = Boundary("3", sys, lambda u,t : u, {2: np.array([0.0, 1.0])}, {2:1})
# n5 = Boundary("4", sys, lambda u,t : u, {3: np.array([1.0, 0.0])}, {3:1})

# test_nodes = [n1, n2, n3, n4, n5]
# #test_nodes = [n1, n2, n3]
# test_nw = Network(test_nodes, test_edges)

