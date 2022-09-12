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
nr_cells = 500
L = 100.0
cfl = 0.9
t = 0
Tfinal = 25.0
cw_order = 3
epsilon = 1.0e-20
alphaPower = 0
CW_type = cppyy.gbl.cweno[m].CW # enum from cweno, i.o.t. move to toml make an if 
RK_type = "ERK_ssp3"
quadOrder = 3

nr_ghosts = math.floor(cw_order/2+1)




cppyy.cppdef(f"""
double  L = {L};
int nr_cells = {nr_cells};
double epsilon = {epsilon};
int alphaPower = {alphaPower};
int cw_order = {cw_order};
double _alpha = {alpha};
double _gamma = {gamma};
""")

cppyy.cppexec(f"""
CLawPsys clawObj(_alpha, _gamma);
Grid G = Grid(0.0, L, nr_cells, {nr_ghosts});
G.setLeftBC(Grid::BC_FREEFLOW);
G.setRightBC(Grid::BC_FREEFLOW);

RecBase<{m}> * Rec_ptr = new cweno<{m}>(G, cw_order, cweno<{m}>::CW, epsilon , alphaPower);
FluxBase<{m}> * numFlux_ptr = new llfFlux<CLawPsys>(clawObj);
numSource<{m}> * Source_ptr = new ZeroSource<{m}>(G);
CLawPsysBCHandler BC;
BCHandler<{m}> * BC_ptr = &BC;
semidiscreteRHS<{m}> * doRHS = new naifRHS<{m}>(G, *Rec_ptr, *numFlux_ptr, *BC_ptr , *Source_ptr);
const int ghosts = doRHS->needsGhosts();
G.setGhosts(ghosts);
ExplicitButcherTableaux * erkTableaux = &getExplicitButcherTableaux<ERKtype>({RK_type});
auto timeIntegrator = explicitRungeKutta<{m}>(*erkTableaux, G ,*doRHS);
""")

# DoF<3, tD> left_U, right_U;
#   for(auto j = 0; j < CLawPsys::M; ++j){
#       left_U[j] = left[j];
#       right_U[j] = right[j];
#   }
#   BC.left_bc =  std::vector<DoF<3, tD>>(cw_order, left_U);
#   BC.right_bc =  std::vector<DoF<3, tD>>(cw_order, right_U);

cppyy.cppdef("""
void set_bc(bool left, std::vector<double> U_std){
  DoF<3, tD> U;
  for(auto j = 0; j < CLawPsys::M; ++j){
      U[j] = U_std[j];
  }
  if (left)
    BC.left_bc =  std::vector<DoF<3, tD>>(cw_order, U);
  else
    BC.right_bc =  std::vector<DoF<3, tD>>(cw_order, U);
}
void set_bc(std::vector<double> left, std::vector<double> right){
  set_bc(true, left);
  set_bc(false, right);
}
""")

physical_range = range(cppyy.gbl.G.beginPhysical(),cppyy.gbl.G.endPhysical()+1 )
dx = cppyy.gbl.G.getDx()
quadrature = gaussRule()
quadrature.chooseOrder(quadOrder)
print(f"Number of nodes in Quadrature is set to: {quadrature.getNNodes()}")

def get_dt(us):
    physical_range
    lambda_max = [cppyy.gbl.clawObj.getLambdaMax(u, 0.0, 0.0) for u in us]
    return cfl*dx/max(lambda_max[physical_range[0]:physical_range[-1]])

# pythonise the DofVec such that ranges can be used!
U = DofVec[m](cppyy.gbl.G.getFullSize())
U1 = DofVec[m](cppyy.gbl.G.getFullSize())

u0 = lambda x : [1.0, 0.0, 0.0] if x < 50.0 else [0.125, 0.0, 0.0]
# set initial

for i in physical_range:
    x = cppyy.gbl.G.getXCenter(i)
    #    cppyy.gbl.clawObj.setU0(cppyy.gbl.CLawPsys.PB_SMOOTHFLAT, x, dx, quadrature, U[i])
    u00 = u0(x)
    U[i][0] = u00[0]    
    
#set boundary
cppyy.gbl.set_bc([1.0, 0.0, 0.0], [0.125, 0.0, 0.0])

fig = plt.figure()
    
# tI = cppyy.gbl.timeIntegrator
# tI.setCFL(cfl)

class cweno3():
    def __init__(self, N, L, cfl, alpha, gamma, init):
        self.N = N
        self.L = L
        self.cfl = cfl
        self.alpha = alpha
        self.gamma = gamma
        self.tI = cppyy.gbl.timeIntegrator
        self.U = DofVec[m](cppyy.gbl.G.getFullSize())
        self.U1 = DofVec[m](cppyy.gbl.G.getFullSize())
        self.physical_range = range(cppyy.gbl.G.beginPhysical(),cppyy.gbl.G.endPhysical()+1 )
        self.tI.setCFL(cfl)
        for i in physical_range:
            x = cppyy.gbl.G.getXCenter(i)
            # cppyy.gbl.clawObj.setU0(cppyy.gbl.CLawPsys.PB_SMOOTHFLAT, x, dx, quadrature, U[i])
            u00 = init(x)
            self.U[i][0] = u00[0]

    def get_dt(self):
        lambda_max = [cppyy.gbl.clawObj.getLambdaMax(u, 0.0, 0.0) for u in self.U]
        return cfl*dx/max(lambda_max[self.physical_range[0]:self.physical_range[-1]])
    def set_dt(self, dt):
        self.tI.setDt(0.1)
    def setDtMax(self, dt):
        self.tI.setDtMax(dt)
    
    def set_both_bc(self, l, r):
        cppyy.gbl.set_bc(l, r)
    def set_bc(self, left_or_right, u):
        cppyy.gbl.set_bc(left_or_right, u)
    # take care the dt has to be set already when calling advance
    def advance(self, t):
        dt = self.tI.advance(self.U, t, self.U1)
        self.t = t + dt
        self.U, self.U1 = self.U1, self.U

# test for the cweno single dynamics
# cw = cweno3(nr_cells, L, cfl, alpha, gamma, u0)
# while (t < Tfinal):
#     cw.setDtMax(Tfinal - t)
#     # tI.setDtMax(Tfinal - t)
#     print(t)
#     cw.set_dt(0.1)
#     cw.set_both_bc([min(1+ t/2, 2), 0.0, 0.0], [min(0.125 + t/2, 1), 0.0, 0.0])
#     cw.advance(t)
#     t = cw.t
#     # dt = tI.advance(U,t, U1)
#     # t = t + dt
#     # U, U1 = U1, U
#     plt.plot([cw.U[i][0] for i in  cw.physical_range])
#     plt.pause(0.05)
#     plt.clf()
# plt.show()



#diamond network


dx = 0.1
u0 = lambda x : [1.0, 0.0, 0.0] if x < 5.0 else [0.5, 0.0, 0.0]
u1 = lambda x : [0.5, 0.0, 0.0]
diamond_edges = {
    (0,1) : Pipe(0, cweno3(0, math.floor(10/dx), 10, cfl, alpha, gamma, u0)),
    (1,2) : Pipe(1, cweno3(1, math.floor(20/dx), np.sqrt(250), cfl, alpha, gamma, u1)),
    (2,3) : Pipe(2, cweno3(2, math.floor(np.sqrt(250)/dx), np.sqrt(250), cfl, alpha, gamma, u1)),
    (3,4) : Pipe(3, cweno3(3, math.floor(10/dx), 10, cfl, alpha, gamma, u1)),
    (5,3) : Pipe(4, cweno3(4, math.floor(20/dx), 20, cfl, alpha, gamma, u1)),
    (1,5) : Pipe(5, cweno3(5, math.floor(20/dx), 20, cfl, alpha, gamma, u1)),
    (2,5) : Pipe(6, cweno3(6, math.floor(3*np.sqrt(200)/2/dx), 3*np.sqrt(200)/2, cfl, alpha, gamma, u1)),
}


sys = Isoeuler_system(1.0, 1.0, 0.0)

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


n1 = Boundary("0", sys,
              {0: np.array([1,0])},
              {0:0},
              coron_db_left)
n2 = Junction("1", sys,
              {0 : np.array([-1,0]),
               1 : np.array([1,0.5])/math.sqrt(1.25),
               5 : np.array([1,-1])/math.sqrt(2)},
              {0:1, 1:0, 5:0})
n3 = Junction("2", sys,
              {1 : np.array([-1,-0.5])/math.sqrt(1.25),
               2 : np.array([1,-0.5])/math.sqrt(1.25),
               6 : np.array([-1,0])},
              {1:1, 2:0, 6:0})
n4 = Junction("3", sys,
              {2 : np.array([-1,0.5])/math.sqrt(1.25),
               3 : np.array([1,0]),
               4 : np.array([-1,-1])/math.sqrt(2)},
              {2:1, 3:0, 4:1})
n5 = Boundary("4", sys,
              {3: np.array([-1,0])},
              {3:1},
              coron_db_right)
n6 = Junction("5", sys,
              {4 : np.array([1,1])/math.sqrt(2),
               5 : np.array([-1,1])/math.sqrt(2),
               6 : np.array([1,0])},
              {4:0, 5:1, 6:1})

diamond_nodes = [n1, n2, n3, n4, n5, n6]
diamond_nw = Network(diamond_nodes, diamond_edges)
# main fctn:
fig = plt.figure()
axs  = fig.add_subplot(111,  projection='3d')
# fig.show()
# fig.canvas.draw()
fig2 = plt.gcf()
#axs21  = fig2.add_subplot(211)
#axs22  = fig2.add_subplot(212)
#fig2.show()

fig.suptitle('Density in a Diamond Network, t = 0', fontsize=20)
plt.xlabel('x')
plt.ylabel('y')

#lines, min_val, max_val = plot_nw(diamond_nw, fig, axs)
#cb = fig.colorbar(lines, ax = axs, orientation='horizontal')
#cb.set_label('Density')


t = 0
Tfinal = 200
plot_nr = 0
i = 10000
while t < Tfinal:
    diamond_nw.advance()
    print(t)
    t = diamond_nw.t
    plot_nr += 1
    if plot_nr == 5:
        #lines, min_val, max_val = plot_nw(diamond_nw, fig, axs)
        plot_nw3d(diamond_nw, fig, axs)
        fig.suptitle(f'Density in a Diamond Network, t = {t:.2f}', fontsize=20)
        plt.pause(0.01)
        plot_nr = 0
        plt.savefig('./img3d/'+str(i)+'.png')
        plt.cla()
        i += 1

# animation('out.avi')


