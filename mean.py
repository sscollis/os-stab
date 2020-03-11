import matplotlib.pyplot as plt
import numpy as np
mean = np.loadtxt("fort.10")
plt.title('Mean flow profile')
plt.ylabel(r"$U(y)$, $U''(y)$")
plt.xlabel(r"$y$")
plt.plot(mean[:,0],mean[:,1],'b-',label=r"$U(y)$")
plt.plot(mean[:,0],mean[:,2],'r-',label=r"$U''(y)$")
plt.legend()
plt.savefig("mean.png")
plt.show(block=False)
