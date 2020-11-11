import matplotlib
import numpy as np

matplotlib.use('Agg')
from matplotlib import pyplot as plt

A = np.loadtxt("lapl.txt", usecols=(5, 9))
B = np.loadtxt("laplv.txt", usecols=(5, 9))

fig = plt.figure(1, frameon=False)
fig.clf()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
ax.plot(np.log2(A[:,1]), A[:,0], marker="o", color="r", ls="-", label="Lapl")
ax.plot(np.log2(B[:,1]), B[:,0], marker="s", color="b", ls="-", label="Lapl-V")
ax.set_xticklabels("$2^{%d}$" % x for x in ax.get_xticks())
ax.legend(loc="upper right", fontsize=9)
ax.set_ylabel(r"$\beta_{FP}$ [Gflop/s]")
ax.set_xlabel(r"$L$")
ax.set_ylim(0, 28)
ax.set_xlim(right=15)
fig.savefig("B.png", transparent=True, dpi=300)
