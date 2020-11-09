import matplotlib
import numpy as np

matplotlib.use('Agg')
from matplotlib import pyplot as plt

A = np.loadtxt("aos.txt", usecols=(5, 13))
B = np.loadtxt("soa-vl16.txt", usecols=(8, 16))
C = np.loadtxt("soa-v.txt", usecols=(8, 16))

fig = plt.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
ax.plot(np.log2(A[:,0]), A[:,1], marker="o", color="r", ls="-", label="AoS")
ax.plot(np.log2(B[:,0]), B[:,1], marker="s", color="b", ls="-", label="SoA, VL=16")
ax.plot(np.log2(C[:,0]), C[:,1], marker="D", color="g", ls="-", label="SoA, intrins.")
ax.set_xticklabels("$2^{%d}$" % x for x in ax.get_xticks())
ax.legend(loc="upper left")
ax.set_ylabel(r"$\beta_{IO}$ [GBytes/s]")
ax.set_xlabel(r"$L$")
fig.savefig("D.png", dpi=300, transparent=True)
