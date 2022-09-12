import numpy as np
from scipy import optimize
from node import *

class Junction(Coupling_node):
    """ Junction coupling node for the isoeuler equations by Colombo, Garavello 2006
            
    Attributes
    ----------
     see the Coupling_node class, no extra attributes are added here

    Methods
    -------
    calc(pipes, t)
      calculates the junction conditions and sets the bounadry conditions of the pipes
    """

    def calc(self, pipes, t):
        """
        Parameters
        ----------
        pipes : dict of pipe_id -> solvers
         dictionary of all pipe solvers connected to this junction node
        t : float
         current time, not used here but is part of the specification, since 
         in general a node may be time dependend
        Returns
        -------
         Nothing
        """
        ids = {}
        bc_vals = []
        i = 0
        for p in pipes:
            pipe_id = p['pipe'].id
            solver = p['pipe'].solver
            ids[pipe_id] = i
            if self.dirs[pipe_id] == 0:
                bc_vals.append(solver.val(0)) 
            else:
                bc_vals.append(solver.val(-2))
                bc_vals[-1][1] *= -1 # settting velocity in the backward direction at outflow, s. below why
            i += 1
        rhos = np.array([u[0] for u in bc_vals])
    
    
        fun = lambda x : self.pde.junction_equations(x, bc_vals, list(self.nus.values()))
        sol = optimize.root(fun, rhos,  method='df-sane',tol=1e-8)
        rhos = sol.x
        # take care for the sign of the velocity! nus are always leaving the node!
        # thus velocity must be taken *(-1) if dir is not zero
        for p in pipes:            
            pipe_id = p['pipe'].id
            if self.dirs[pipe_id] == 0:
                p['pipe'].solver.set_bc(True, self.pde.lax_curve(2, -1, rhos[ids[pipe_id]], bc_vals[ids[pipe_id]] ))
            else:
                bc_vals[ids[pipe_id]][1] *= -1
                p['pipe'].solver.set_bc(False, self.pde.lax_curve(1, 1, rhos[ids[pipe_id]], bc_vals[ids[pipe_id]] ))


            

class Boundary(Coupling_node):
    """ Boundary leaf for the isoeuler equations
            
    Attributes
    ----------
    see the Coupling_node class and in addition:
    
    f(u, t)=id : function
     boundary function, where u is the FV solution of adjacent pipe and t current time
    Methods
    -------
    calc(pipes, t)
      calculates the boundary conditions and sets them in the pipes solver
    """    
    def __init__(self, id, pde, nus, dirs, f=id):
        """
        Parameters
        ----------

        id : int/str
         Unique id of the coupling node
        pde : class
         pde system class, see pde.py for the specifications of
         the required attributes and methods of the pde classes and isoeuler.py
         for an example FV implementation
        nus : dict of pipe_id -> np.array each with length 2
         direction of each pipe relative to the node with its 2-norm 
         encodes the thickness of the adjacent pipe.
        dirs : dict of pipe_id -> int
         indication which end of the pipe is connected to the node, if pipe_id:0
         the pipe with id pipe_id is connected with its left-hand side boundary
         to the node and setting pipe_id:1 indicates that the pipe is connected
         with the right-hand side boundary.
        f(u, t)=id : function
         boundary function, where u is the FV solution of adjacent pipe and t current time
        """
        # todo: check sizes nus and dir should be exactly one
        self.id = id
        self.pde = pde
        self.f = f
        self.nus = nus
        self.dirs = dirs
    def calc(self, pipes_list, t):
        """
        Parameters
        ----------
        pipes : dict of pipe_id -> solvers
         dictionary of all pipe solvers connected to this junction node
        t : float
         current time, used here in the function self.f
        Returns
        -------
         Nothing
        """
        pipe_id = list(self.dirs.keys())[0]
        solver = pipes_list[0]['pipe'].solver
        pipe_id = pipes_list[0]['pipe'].id
        if(self.dirs[pipe_id] == 0):
            solver.set_bc(True, self.f(solver.val(0), t))
        else:
            solver.set_bc(False, self.f(solver.val(-1), t))



def coron_db_left(u, t):
    nu = 0.1
    R = sys.riemann_invariants(u)
    Rminus = R[0]
    Rplus = nu*R[0]
    return sys.conservative(np.array([Rminus, Rplus]))
    
def coron_db_right(u, t):
    nu = 0.1
    R = sys.riemann_invariants(u)
    Rminus = nu*R[1]
    Rplus = R[1]
    return sys.conservative(np.array([Rminus, Rplus]))
