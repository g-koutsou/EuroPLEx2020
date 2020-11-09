import matplotlib
import numpy as np

matplotlib.use('Agg')
from matplotlib import pyplot as plt

A = np.loadtxt("aos.txt", usecols=(5, 13))

fig = plt.figure(1, frameon=False)
fig.clf()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
ax.plot(np.log2(A[:,0]), A[:,1], marker="o", color="r", ls="-", label="AoS")
ax.set_xticklabels("$2^{%d}$" % x for x in ax.get_xticks())
ax.legend(loc="upper left")
ax.set_ylabel(r"$\beta_{IO}$ [GBytes/s]")
ax.set_xlabel(r"$L$")
fig.savefig("A.png", transparent=True, dpi=300)
