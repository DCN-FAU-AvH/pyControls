import math
import cppyy, cppyy.ll
cppyy.add_include_path('../lib/claw1darena/src')
cppyy.add_include_path('../lib/claw1darena/build')
cppyy.include('claw1darena_all.hpp')

from cppyy.gbl import CLawPsys, Grid, cweno, RecBase, llfFlux, CLawSWEBCHandler, ZeroSource, naifRHS, Euler, DofVec, gaussRule
cppyy.load_library('libboost_program_options')
cppyy.load_library('../lib/claw1darena/build/libshared_claw1darena')

import matplotlib.pyplot as plt

# move to toml settings
m = 3
alpha = 1.0
gamma = 1.0
w = 0.1337
nr_cells = 500
L = 100.0
cfl = 0.9
t = 0
Tfinal = 40.0
cw_order = 3
epsilon = 1.0e-20
alphaPower = 0
CW_type = cppyy.gbl.cweno[m].CW # enum from cweno, i.o.t. move to toml make an if 
RK_type = "ERK_ssp3"
quadOrder = 3

nr_ghosts = math.floor(cw_order/2+1)



def setup_cweno_psys(id, alpha, gamma, L, nr_cells):
    cppyy.cppexec(f"""
    auto cweno3{id} = cweno3_psys({alpha}, {gamma}, {w}, {L}, {nr_cells}, {nr_ghosts}, {epsilon}, {alphaPower} )
    """)


# todo add destructor to cweno3_psys

# the id argument generates expressions for all c++ objects like G1, G2, etc because here
# everything is global namespace
class cweno3():
    def __init__(self, id, N, L, cfl, alpha, gamma, init):
        self.id = id
        self.N = N
        self.L = L
        self.t = 0
        self.cfl = cfl
        self.alpha = alpha
        self.gamma = gamma
        setup_cweno_psys(id, alpha, gamma, L, N)
        self.tI = eval(f"cppyy.gbl.cweno3{id}.get_timeIntegrator()")
        self.G = eval(f"cppyy.gbl.cweno3{id}.get_Grid()")
        self.clawObj = eval(f"cppyy.gbl.cweno3{id}.clawObj")        
        self.U = DofVec[m](self.G.getFullSize())
        self.U1 = DofVec[m](self.G.getFullSize())
        self.physical_range = range(self.G.beginPhysical(), self.G.endPhysical()+1 )
        self.tI.setCFL(cfl)
        #todo separate function for initilisation and for all variables
        for i in self.physical_range:
            x = self.G.getXCenter(i)
            # cppyy.gbl.clawObj.setU0(cppyy.gbl.CLawPsys.PB_SMOOTHFLAT, x, dx, quadrature, U[i])            
            u00 = init(x)
            self.U[i][0] = u00[0]
    def val(self, i):
        res = []
        for j in range(0,2): 
            res.append(self.U[self.physical_range[i]][j])
        return res
    def get_dt(self):
        lambda_max = [self.clawObj.getLambdaMax(u, 0.0, 0.0) for u in self.U]
        dx = self.G.getDx()
        return cfl*dx/max(lambda_max[self.physical_range[0]:self.physical_range[-1]])
    def set_dt(self, dt):
        self.tI.setDt(dt)
    def setDtMax(self, dt):
        self.tI.setDtMax(dt)
    
    def set_both_bc(self, l, r):
        eval(f"cppyy.gbl.cweno3{self.id}.set_both_bc(l, r)")        
    def set_bc(self, left_or_right, u):
        eval(f"cppyy.gbl.cweno3{self.id}.set_bc(left_or_right, u)")        

    # take care the dt has to be set already when calling advance
    def advance(self, t):
        dt = self.tI.advance(self.U, t, self.U1)
        self.t = t + dt
        self.U, self.U1 = self.U1, self.U

# u0 = lambda x : [1.0, 0.0, 0.0] if x < 50.0 else [0.125, 0.0, 0.0]
# cw = cweno3(0, nr_cells, L, cfl, alpha, gamma, u0)
# while (t < Tfinal):
#     with cppyy.ll.signals_as_exception():
#         cw.setDtMax(Tfinal - t)
#         print(t)
#         cw.set_both_bc([1.5, 0.0, 0.0], [0.5, 0.0, 0.0])
#         cw.advance(t)
#         t = cw.t
#         t = t + 0.1#dt
#         # U, U1 = U1, U
#         plt.plot([cw.U[i][0] for i in  cw.physical_range])
#         plt.pause(0.1)
#         plt.clf()
# plt.show()




