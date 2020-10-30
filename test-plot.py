import numpy as np
from matplotlib import pyplot as plt

x = np.linspace(-np.pi, np.pi, 25)
y = np.tanh(x)

fig = plt.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
ax.plot(x, y, ls="", marker="s")
ax.set_xlabel(r"$\theta$")
ax.set_ylabel(r"tanh($\theta$)")
fig.savefig("test-plot.pdf")
