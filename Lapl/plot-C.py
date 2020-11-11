import matplotlib
import numpy as np

matplotlib.use('Agg')
from matplotlib import pyplot as plt

A = np.loadtxt("lapl.txt", usecols=(5, 9))
B = dict()
for vl in [4, 8, 16, 32, 64]:
    B[vl] = np.loadtxt("laplv-VL%d.txt" % vl, usecols=(5, 9))

fig = plt.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
### ax.set_yscale("log")
ax.plot(np.log2(A[:,1]), A[:,0], marker="o", color="r", ls="-", label="Lapl")
for vl in sorted(B):
    ax.plot(np.log2(B[vl][:,1]), B[vl][:,0], marker="s", ls="-", label="Lapl-V, VL=%d" % vl)
ax.set_xticklabels("$2^{%d}$" % x for x in ax.get_xticks())
ax.legend(loc="upper right", fontsize=9)
ax.set_ylabel(r"$\beta_{FP}$ [Gflop/s]")
ax.set_xlabel(r"$L$")
ax.set_ylim(0, 28)
ax.set_xlim(right=15)
fig.savefig("C.png", dpi=300, transparent=True)
