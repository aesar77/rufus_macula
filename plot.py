import matplotlib
matplotlib.rcParams['font.family'] = 'serif'
import matplotlib.pyplot as plt
import numpy as np

dt = 0.01
for n in range(399):

    t = n * dt
    print(t)
    data = np.loadtxt(f"output3d_rho_{n}.csv", delimiter=",")
    # print(data[0])
    fig, ax = plt.subplots()
    cax =  ax.imshow(data, cmap="inferno", extent=[0,1,0,1])
    ax.set_title(fr'$\rho$, t={t:.2f}')

# Add colorbar, make sure to specify tick locations to match desired ticklabels
    cbar = fig.colorbar(cax)
    # cbar.ax.set_yticklabels(['< -0.1', '0', '> 0.1']) 
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(f"rho3d_{n:03d}.png", bbox_inches="tight")
    plt.close()

