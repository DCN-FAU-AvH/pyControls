import math
import numpy as np
import matplotlib.pyplot as plt
import network as nw

# import iseuler stuff
from isoeuler import llf, pde, nodes
# diamond network

sys = pde.Isoeuler_system(1.0, 1.0, 0.0)

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

dx = 0.5
u0 = lambda x : [1.0, 0.0] if x < 5.0 else [0.5, 0.0]
u1 = lambda x : [0.5, 0.0]

cfl = 0.5
diamond_edges = {
    (0,1) : nw.Pipe(0, llf.iso_llf(0, math.floor(10/dx), 10, cfl, sys, u0)),
    (1,2) : nw.Pipe(1, llf.iso_llf(1, math.floor(20/dx), np.sqrt(250), cfl, sys, u1)),
    (2,3) : nw.Pipe(2, llf.iso_llf(2, math.floor(np.sqrt(250)/dx), np.sqrt(250), cfl, sys, u1)),
    (3,4) : nw.Pipe(3, llf.iso_llf(3, math.floor(10/dx), 10, cfl, sys, u1)),
    (5,3) : nw.Pipe(4, llf.iso_llf(4, math.floor(20/dx), 20, cfl, sys, u1)),
    (1,5) : nw.Pipe(5, llf.iso_llf(5, math.floor(20/dx), 20, cfl, sys, u1)),
    (2,5) : nw.Pipe(6, llf.iso_llf(6, math.floor(3*np.sqrt(200)/2/dx), 3*np.sqrt(200)/2, cfl, sys, u1)),
}

n1 = nodes.Boundary("0", sys,
                 {0: np.array([1,0])},
                 {0:0},
                 coron_db_left)
n2 = nodes.Junction("1", sys,
                 {0 : np.array([-1,0]),
                  1 : np.array([1,0.5])/math.sqrt(1.25),
                  5 : np.array([1,-1])/math.sqrt(2)},
                 {0:1, 1:0, 5:0})
n3 = nodes.Junction("2", sys,
                 {1 : np.array([-1,-0.5])/math.sqrt(1.25),
                  2 : np.array([1,-0.5])/math.sqrt(1.25),
                  6 : np.array([-1,0])},
                 {1:1, 2:0, 6:0})
n4 = nodes.Junction("3", sys,
                 {2 : np.array([-1,0.5])/math.sqrt(1.25),
                  3 : np.array([1,0]),
                  4 : np.array([-1,-1])/math.sqrt(2)},
                 {2:1, 3:0, 4:1})
n5 = nodes.Boundary("4", sys,
                 {3: np.array([-1,0])},
                 {3:1},
                 coron_db_right)
n6 = nodes.Junction("5", sys,
                 {4 : np.array([1,1])/math.sqrt(2),
                  5 : np.array([-1,1])/math.sqrt(2),
                  6 : np.array([1,0])},
                 {4:0, 5:1, 6:1})

diamond_nodes = [n1, n2, n3, n4, n5, n6]

diamond_nw = nw.Network(diamond_nodes, diamond_edges)
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


# nw = Network(nodes_map, dynamics_map )
    

# animation('out.avi')

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
        nw.plot_nw3d(diamond_nw, fig, axs)
        fig.suptitle(f'Density in a Diamond Network, t = {t:.2f}', fontsize=20)
        plt.pause(0.01)
        plot_nr = 0
        plt.savefig('./img3d/'+str(i)+'.png')
        plt.cla()
        i += 1
#    fig.colorbar(line, ax=axs)
#    fig.canvas.draw()
    # print(t)
    # plt.cla() 
    # axs21.plot([test_edges[(0,1)].solver.U[i][0] for i in  test_edges[(0,1)].solver.physical_range])
    # axs22.plot([test_edges[(0,2)].solver.U[i][0] for i in  test_edges[(0,1)].solver.physical_range])
    #plt.pause(0.1)
