from gaslib import *

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
#    fig.colorbar(line, ax=axs)
#    fig.canvas.draw()
    # print(t)
    # plt.cla() 
    # axs21.plot([test_edges[(0,1)].solver.U[i][0] for i in  test_edges[(0,1)].solver.physical_range])
    # axs22.plot([test_edges[(0,2)].solver.U[i][0] for i in  test_edges[(0,1)].solver.physical_range])
    #plt.pause(0.1)
