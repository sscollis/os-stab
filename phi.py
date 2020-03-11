import matplotlib.pyplot as plt
import numpy as np
phi = np.loadtxt("fort.11")
plt.title('Eigenfunction')
plt.ylabel(r"$\Re(\phi(y))$, $\Im(\phi(y))$")
plt.xlabel(r"$y$")
plt.plot(phi[:,0],phi[:,1],'b-',label=r"$\Re(\phi(y))$")
plt.plot(phi[:,0],phi[:,2],'r-',label=r"$\Im(\phi(y))$")
plt.legend()
plt.savefig("phi.png")
plt.show(block=False)
