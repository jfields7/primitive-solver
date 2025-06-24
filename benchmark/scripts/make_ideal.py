import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mplot

plt.style.use("plots.mplstyle")

# Data information
fname = "../bin/ideal.dat"
nn = 100
nT = 100

# Conversion factors from code units to CGS
conv_temp = 10903137627740.92
conv_dens = 6.176269145886162e+17

# Load the data
data = np.genfromtxt(fname, usecols=(0,1,3))
rho = data[:,0]
T = data[:,1]
error = data[:,2]
iters = np.genfromtxt(fname, usecols=(2))
success = np.genfromtxt(fname, usecols=(4), dtype=np.uint32)

mask = success == 0

iters[mask] = np.nan
error[mask] = np.nan

# Reshape the data
rho = np.reshape(rho, newshape=(nT, nn))*conv_dens
T = np.reshape(T, newshape=(nT, nn))*conv_temp
error = np.reshape(error, newshape=(nT, nn))
iters = np.reshape(iters, newshape=(nT, nn))
success = np.reshape(iters, newshape=(nT, nn))

print(np.any(np.logical_not(success)))

#iters[np.logical_not(success)] = 100

colors1 = plt.cm.gist_rainbow_r(np.linspace(0.05, 0.2, 600))
colors2 = plt.cm.gist_rainbow_r(np.linspace(0.2, 0.75, 500))
colors3 = plt.cm.gist_rainbow_r(np.linspace(0.75, 1, 1000))
cmap = colors.LinearSegmentedColormap.from_list('colormap',np.vstack((colors1, colors2, colors3)))

#omap = mplot.colormaps['gist_rainbow']
#rmap = omap.reversed()

cs = plt.pcolormesh(rho, T, iters, shading='gouraud',
                    norm=colors.Normalize(vmin=1, vmax=30), cmap=cmap)
#cs.cmap.set_under('w')
cs.cmap.set_over('k')
plt.colorbar(extend='max', label='# Iterations')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\rho~[\mathrm{g/cm}^3]$')
plt.ylabel(r'$T~[\mathrm{K}]$')

plt.savefig('ideal_iter.png')
plt.savefig('ideal_iter.pdf')

plt.show()

mask = error < 1e-15
error[mask] = 1e-15

# Get the average number of iterations and the average error
avg_iters = np.nanmean(iters)
avg_err = 10.0**(np.nanmean(np.log10(error)))

print(f"Avg. iters: {avg_iters}")
print(f"Avg. err: {avg_err}")


colors1 = plt.cm.viridis(np.linspace(0.1, 1.0, 800))
colors2 = plt.cm.inferno_r(np.linspace(0.067, 0.9, 800))
cmap = colors.LinearSegmentedColormap.from_list('colormap',np.vstack((colors1, colors2)))

cs = plt.pcolormesh(rho, T, error, shading='gouraud',
                    norm=colors.LogNorm(vmin=1e-15, vmax=5e-4), cmap=cmap)
cs.cmap.set_over('w')
plt.colorbar(extend='max', label='Relative error')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\rho~[\mathrm{g/cm}^3]$')
plt.ylabel(r'$T~[\mathrm{K}]$')

plt.savefig('ideal_err.png')
plt.savefig('ideal_err.pdf')

plt.show()
