import numpy as np
from itertools import chain
import matplotlib.pyplot as plt

class fv():
    """
    Pipe dynamics via finite volumes base class

    Attributes
    ----------
    id : int/str
     id of the pipe dynamics object
    N : int
     number of cells
    L : float
     length of the pipe
    t : float
     current time
    ghosts : int
     number of required ghost cells
    cfl : float
     CFL number
    dt_max : float
     maximal available dt (default = Inf), required for e.g. cweno from claw1darena
    pde : pde
     pde object, cf. pde.py or isoeuler for an example
    U : np.array
     fv solution data including ghost cells
    U1 : np.array
     fv solution data including ghost cells temporary, advance method sets U1 to the updated values and then swaps U and U1
    physical_range : range
     indices of the solution U (and U1) without ghost cells
    dx : float
     cell width of the uniform grid
    dt : float
     current time step used in advance

    Methods
    -------
    val(i)
     returns the value of the solution U of the physical grid, i.e. val(i) evaluates U(ghosts + i)
    get_dt() 
     returns the timestep according to the cfl condition
    set_dt(dt)
     sets the timestep, i.e. sets self.dt
    setDtMax(dt)
     sets the maximal timestep, i.e. sets self.max_dt
    set_both_bc(l, r)
     sets the lhs ghosts constantly to l and the rhs to r
    set_bc(left_or_right, u)
     if left_or_right = True, sets the left ghosts to u and the rhs otherwise
    nf()
     evaluates the numerical flux function, *has to be implemented in child class!*, just raises an exception here
    source()
     evaluates the numerical source discretisation, *has to be implemented in child class!*, just raises an exception here
    advance(t)
     updates the solution by performing one timestep using the first order Euler step and numerical flux discretisation,
     dt has to be set before calling advance
    """

    def __init__(self, id, N, L, cfl, pde, init, ghosts=1):
        """
        Parameters
        ----------
        id : int/str
         id of the pipe dynamics object
        N : int
         number of cells
        L : float
         length of the pipe
        cfl : float
         CFL number
        pde : pde
         pde system object, see pde.py and isoeuler.py for an example
        init : function
         initial conditions
        ghosts : int
         number of required ghost cells
        """
        self.id = id
        self.N = N
        self.L = L
        self.t = 0
        self.ghosts = ghosts
        self.cfl = cfl
        self.dt_max = float('Inf')
        self.pde = pde
        self.U = np.zeros((pde.M, N + 2*ghosts))
        self.U1 = np.zeros((pde.M, N + 2*ghosts))
        self.physical_range = range(ghosts, N+ghosts)
        self.dx = float(L)/N
        for i in self.physical_range:
            self.U[:, i] = init(self.dx/2 + i*self.dx)
        self.set_both_bc(self.U[:, ghosts], self.U[:, N])    
        self.dt = self.get_dt()
        
    def val(self, i):
        """returns the value of the solution U of the physical grid, i.e. val(i) evaluates U(ghosts + i)

        Parameters
        ----------
        i : int

        Returns
        _______
        
        """
        return self.U[:, self.ghosts + i]
    
    def get_dt(self):
        """returns the timestep size according to the cfl condition

        Returns
        _______
        float
         timestep size according to the cfl condition
        """
        max_ev = 0.0
        for i in self.physical_range:
            l = abs(self.pde.ev(self.U[:, i]))
            max_ev = max(max_ev, l.max())

        return self.cfl*self.dx/max_ev
    
    def set_dt(self, dt):
        """sets the timestep, i.e. sets self.dt

        Parameters
        ----------
        dt : float
         timestep size
        """
        self.dt = dt
        
    def setDtMax(self, dt):
        """sets the maximal timestep, i.e. sets self.max_dt

        Parameters
        ----------
        dt : float
         maximal timestep size
        """
        self.dt_max = dt

    def set_both_bc(self, l, r):
        """sets the lhs ghosts constantly to l and the rhs to r

        Parameters
        ----------
        l : np.array
         of length pde.M, lhs bc ghost cell values
        r : np.array
         of length pde.M, rhs bc ghost cell values
        """
        for i in range(0, self.ghosts):
            self.U[:,i] = l
        for i in range(self.N+self.ghosts,self.N+2*self.ghosts):
            self.U[:,i] = r
            
    def set_bc(self, left_or_right, u):
        """if left_or_right = True, sets the left ghosts to u and the rhs otherwise
        
        Parameters
        ----------
        left_or_right : bool
         if True, sets the left cell, otherwise right cell
        u : np.array
         of length pde.M, bc ghost cell values
        """
        if left_or_right:
            for i in range(0, self.ghosts):
                self.U[:,i] = u
        else:
            for i in range(self.N+self.ghosts,self.N+2*self.ghosts):
                self.U[:,i] = u

    def nf(self):
        """evaluates the numerical flux function for self.U, *has to be implemented in child class!*, just raises an exception here              

        Returns
        _______
        np.array
         array of numerical fluxes nfs, nfs[i] is the lhs numerical flux of the physical(!) cell i and nfs[i+1] the rhs        
        """
        raise RuntimeError("nf has to be implemented (in a child class)")

    def source(self):
        """evaluates the numerical source discretisation for self.U, *has to be implemented in child class!*, just raises an exception here

        Returns
        _______
        np.array
         array of disc. sources, that is indexed as the physical_range (!)
        """
        raise RuntimeError("source has to be implemented (in a child class)")

    def advance(self, t):
        """updates the solution by performing one timestep using the first order Euler step and numerical flux discretisation,
        dt has to be set before calling advance

        Parameters
        ----------
        t : float
         current time

        Returns
        _______
        float
         current dt        
        """
        nfs = self.nf()
        src = self.source()
        for i in self.physical_range:
            self.U1[:, i] = self.U[:, i] - self.dt/self.dx*(nfs[:,i+1] - nfs[:,i]) + src[:, i]*self.dt
        self.U1[:, range(0, self.ghosts)] = self.U[:, range(0, self.ghosts)]
        self.U1[:, range(self.N+self.ghosts,self.N+2*self.ghosts)] = \
            self.U[:, range(self.N+self.ghosts,self.N+2*self.ghosts)]
        self.t = t + self.dt
        self.U, self.U1 = self.U1, self.U
        return self.dt

