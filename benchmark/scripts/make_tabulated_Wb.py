import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mplot

plt.style.use("plots.mplstyle")

# Data information
#fname = "../bin/dd2_Wb.dat"
fname = "../bin/ls220_wb.dat"
nbeta = 100
nW = 100

# Load the data
data = np.genfromtxt(fname, usecols=(0,1,3))
W = data[:,0]
ibeta = data[:,1]
error = data[:,2]
iters = np.genfromtxt(fname, usecols=(2))
success = np.genfromtxt(fname, usecols=(4), dtype=np.uint32)

print(np.sum(success)/success.shape[0])

mask = success == 0

iters[mask] = np.nan
error[mask] = np.nan

# Reshape the data
W = np.reshape(W, newshape=(nbeta, nW)) - 1.0
ibeta = np.reshape(ibeta, newshape=(nbeta, nW))
error = np.reshape(error, newshape=(nbeta, nW))
iters = np.reshape(iters, newshape=(nbeta, nW))
success = np.reshape(success, newshape=(nbeta, nW))

#print(np.any(np.logical_not(success)))

#iters[np.logical_not(success)] = 100

colors1 = plt.cm.gist_rainbow_r(np.linspace(0.05, 0.2, 600))
colors2 = plt.cm.gist_rainbow_r(np.linspace(0.2, 0.75, 500))
colors3 = plt.cm.gist_rainbow_r(np.linspace(0.75, 1, 1000))
cmap = colors.LinearSegmentedColormap.from_list('colormap',np.vstack((colors1, colors2, colors3)))

#omap = mplot.colormaps['gist_rainbow']
#rmap = omap.reversed()

mask = success == 1

cs = plt.pcolormesh(W, ibeta, iters, shading='gouraud',
                    norm=colors.Normalize(vmin=1, vmax=30), cmap=cmap)
#cs.cmap.set_under('w')
cs.cmap.set_over('k')
plt.colorbar(extend='max', label='# Iterations')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1e-3, 1e3])
plt.ylim([1e-10, 1e10])
plt.xlabel(r'$W-1$')
plt.ylabel(r'$P_\mathrm{mag}/P$')

plt.savefig('ls220_wb_iter.pdf')
plt.savefig('ls220_wb_iter.png')

plt.show()

mask = error < 1e-15
error[mask] = 1e-15
mask = success == 1

colors1 = plt.cm.viridis(np.linspace(0.1, 1.0, 800))
colors2 = plt.cm.inferno_r(np.linspace(0.067, 0.9, 800))
cmap = colors.LinearSegmentedColormap.from_list('colormap',np.vstack((colors1, colors2)))

cs = plt.pcolormesh(W, ibeta, error, shading='gouraud',
                    norm=colors.LogNorm(vmin=1e-15, vmax=5e-4), cmap=cmap)
cs.cmap.set_over('w')
plt.colorbar(extend='max', label='Relative error')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1e-3, 1e3])
plt.ylim([1e-10, 1e10])
plt.xlabel(r'$W-1$')
plt.ylabel(r'$P_\mathrm{mag}/P$')

plt.savefig('ls220_wb_err.pdf')
plt.savefig('ls220_wb_err.png')

plt.show()
