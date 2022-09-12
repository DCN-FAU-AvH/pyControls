from num import *
class iso_llf(fv):
    """
    Pipe dynamics for the isoeuler equations with local Lax Friedrich flux

    Attributes
    ----------
    see parent class fv

    Methods
    -------
    see parent class fv, nf and source are implemented here
    """
    def nf(self):
        """evaluates the numerical flux function for the isoeuler equations using the local Lax Friedrich scheme for self.U

        Returns
        _______
        np.array
         array of numerical fluxes nfs, nfs[i] is the lhs numerical flux of the physical(!) cell i and nfs[i+1] the rhs        
        """
        fluxes = self.pde.flux(self.U)

        nfs = np.zeros(self.U.shape)
        for i in range(self.ghosts, self.N+self.ghosts+1):
            alpha = max(abs(self.pde.ev(self.U[:,i])[0]),
                        abs(self.pde.ev(self.U[:,i])[1]),
                        abs(self.pde.ev(self.U[:,i-1])[0]),
                        abs(self.pde.ev(self.U[:,i-1])[1]))
            nfs[:, i] =  0.5*(fluxes[:,i] + fluxes[:,i-1]) - 0.5*alpha*(self.U[:,i] - self.U[:,i-1])
        return nfs

    def source(self):
        """evaluates the numerical source discretisation for self.U for the isoueuler equations
                
        Returns
        _______
        np.array
         array of disc. sources, that is indexed as the physical_range (!)
        """

        sources = np.zeros(self.U.shape)
        for i in self.physical_range:
            x = self.dx/2 + i*self.dx
            sources[1, i] = -0.5 * (self.pde.theta + self.pde.w(x))*\
                self.U[1, i]*abs(self.U[1, i])/self.U[0, i];
        return sources

# test
# import isoeuler
# sys = isoeuler.Isoeuler_system(1.0, 1.0, 0.0)
# fvol = iso_llf(0, 300, 10, 0.99, sys, lambda x :  [1, 0.0] if x > 5.0 else [0.3, 0.0], 2)
# fvol.set_both_bc([0.3, 0.0], [1.0, 0.0])

# Tfinal = 0.5
# t = 0.0

# while (t < Tfinal):
#         print(t)
#         fvol.dt = fvol.get_dt()
#         dt = fvol.advance(t)
#         t = t + dt
#         plt.plot([fvol.U[0, i] for i in  fvol.physical_range])
#         plt.pause(0.1)
#         plt.clf()
# plt.show()
