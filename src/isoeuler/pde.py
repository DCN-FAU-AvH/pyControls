import math
import numpy as np

class Isoeuler_system(pde):
    """
    The 1D isoeuler pde system u_t + f(u)_x = S(u, x)

    with u = (\rho, \rho v)^T and the pressure law p(\rho) = alpha*rho^gamma
    and friction source term S(u, x) = (theta + w(x))\rhov|\rhov| / \rho

    Attributes
    ----------
    alpha : float
     pressure law constant
    gamma : float
     adiabatic exponent
    theta : float
     friction coefficient
    w(x) : function
     perturbation of the friction coefficient
    M : int
     number of equations in the system, here M=2

    Methods
    -------
    pressure(rho)
     speakes of itself
    flux(u)
     return the flux funciton f
    ev(u)
     returns the eigenvalues of the Jacobian of f
    lax_curve(ith_curve, direction, rho, u)
     return the conservative value along a lax curve
    junction_equation(rhos, u0, nus)
     returns the rhs of the junction equations by Colombo, Garavello 2006
    riemann_invariants(u)
     converts the conservative variables vector u to Riemann invariants vector
    conservative(R)
     converts the Riemann invariants vector R to conservative variables vector
    """
    
    def __init__(self, alpha, gamma, theta=0, w=lambda x : 0.0):
        """
        Parameters
        ----------

        alpha : float
         pressure law constant
        gamma : float
         adiabatic exponent
        theta=0 : float
         friction coefficient
        w(x)=lambda x : 0.0 : function
         perturbation of the friction coefficient

        """
        self.alpha = alpha
        self.gamma = gamma
        self.theta = theta
        self.w = w
        self.M = 2

    def pressure(self, rho):
        """
        Pressure gamma law

        Parameters
        ----------       
        rho : float
         density
        
        Returns
        -------
        float
         pressure p
        """
        return self.alpha*rho**self.gamma

    def flux(self,u):
        """
        The flux funciton f

        Parameters
        ----------       
        u : np.array 
         of length M, a vector of conservative variables
                 
        Returns
        -------
        np.array 
         flux function f
        """
        return np.array([u[1], u[1]*u[1]/u[0] + self.pressure(u[0]) ])

    def ev(self, u):
        """
        The eigenvalues of the Jacobian of f

        Parameters
        ----------       
        u : np.array 
         of length M, a vector of conservative variables
                 
        Returns
        -------
        np.array 
         vector of eigenvalues in ascending order
        """
        v = u[1]/u[0]
        c = math.sqrt(self.gamma * self.pressure(u[0])/u[0])
        return np.array([v - c, v + c])

    def  lax_curve(self, ith_curve, direction, rho, u):
        """
        Return the conservative value along a lax curve

        Parameters
        ----------       
        ith_curve : int
         nr of the lax curve <= M
        direction : int
         of the lax curve, either positive or negative value
        rho : float
         density parameter of the Lax curve
        u : np.array 
         of length M, a vector of conservative variables, foot point of the lax curve

        Raises
        ------
        todo: together with the exceptions in the implementation
                 
        Returns
        -------
        np.array 
         conservative variables vector along the lax curve ~ith_curve~ in direction ~direction~
         with parameter rho, i.e. the resulting value has density rho
        """
        int_c_over_rho = 0.0
        lax_u = np.array([0.0, 0.0])
        p = lambda rho : self.pressure(rho)
  

        if self.gamma == 1.0:
            int_c_over_rho = rho * math.sqrt(self.alpha) * math.log(rho / u[0])
        elif self.gamma > 1.0:
            int_c_over_rho = rho * math.sqrt(self.alpha * self.gamma) * 2.0 / (self.gamma - 1.0) * \
                (rho ** (0.5 * (self.gamma - 1.0)) - u[0] ** (0.5 * (self.gamma - 1.0)))
        else:
            # todo: exception
            pass
    
        sqrt1 = math.sqrt(rho / u[0] * abs(rho - u[0]) * abs(p(rho) - p(u[0])))

        lax_u[0] = rho
        #first term in the resulting momentum. The second comes right after
        lax_u[1] = rho * u[1] / u[0]
        if direction > 0:
            #forward lax curves
            if ith_curve == 1:
                if rho < u[0]:
                    lax_u[1] = lax_u[1] - int_c_over_rho
                else:
                    lax_u[1] = lax_u[1] -  sqrt1
                    
            elif ith_curve == 2:
                if rho < u[0]:
                    lax_u[1] = lax_u[1] - sqrt1
                else:
                    lax_u[1] = lax_u[1] + int_c_over_rho
            else:
                #throw(DomainError(ith_curve, "Error: ith_curve not defined in for the psystem"))
                pass
        elif direction < 0:
            #backward lax curve
            if ith_curve == 1:
                if rho < u[0]:
                    lax_u[1] = lax_u[1] + sqrt1
                else:
                    lax_u[1] = lax_u[1] - int_c_over_rho

            elif ith_curve == 2:
                if rho < u[0]:
                    lax_u[1] = lax_u[1] + int_c_over_rho
                else:
                    lax_u[1] = lax_u[1] + sqrt1
        
            else:
                #throw(DomainError(ith_curve, "Error: ith_curve not defined in for the psystem"))
                pass
        else:
            #throw(DomainError(direction, "Error: Direction not defined in for the psystem"))
            pass
        return lax_u
   

    def junction_equations(self, rhos, u_0, nus):
        """
        The rhs of the junction equations by Colombo, Garavello 2006

        Parameters
        ----------       
        rhos : list
         of same legnth as nus, a vector of denisities
        u_0 : list
         list of foot_point conservative value vectors
        nus : list
         of directional vectors of the junction node
                 
        Returns
        -------
        np.array 
         rhs of the junction equations by Colombo, Garavello 2006
        """
        if (rhos <= 0).any():
            return float('inf')
        # todo: chekc that u_0, rhos and nus are same length!
        nus_norms = [np.linalg.norm(nu) for nu in nus]
        K = len(nus)
        qs =  np.zeros(K)
        ps =  np.zeros(K)

        for i in range(0,K):
            lc = self.lax_curve(2, -1, rhos[i], u_0[i])
            qs[i] = lc[1]
            ps[i] = self.pressure(lc[0]) + np.power(qs[i], 2)/lc[0]

        result = np.zeros(K)
        result[0] = np.dot(nus_norms, qs)
        for i in range(1,K):
            result[i] = ps[i-1] - ps[i]
        return result

    def riemann_invariants(self, u):
        """
        Converts the conservative variables vector u to Riemann invariants vector

        Parameters
        ----------       
        u : np.array 
         of length M, a vector of conservative variables
                 
        Returns
        -------
        np.array 
         Riemann invariants
        """
        p = self.pressure(u[0])
        return np.array( [math.log(p) - math.sqrt(self.alpha)*u[1]/p,
                              math.log(p) + math.sqrt(self.alpha)*u[1]/p])

    def conservative(self, R):
        """
        Converts the Riemann invariants vector R to conservative variables vector

        Parameters
        ----------       
        u : np.array 
         of length M, a vector of Riemann invariants
                 
        Returns
        -------
        np.array 
         Conservative variables
        """
        p = math.exp(0.5*(R[0] + R[1]))
        rho = p/self.alpha; # todo correct that according to paper when using compr factor
        q = p*(R[1] - R[0])/2/math.sqrt(self.alpha)
        return np.array([rho, q])


