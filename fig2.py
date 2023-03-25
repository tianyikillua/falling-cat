import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import degree
from FallingCat import FallingCat

JI = 0.25
alpha = 30*degree

plt.figure(figsize=(5,7))

c = FallingCat(JI, alpha)
t = c.theta/degree
psi = c.lean()/degree
gamma = c.bend()/degree
phi = c.twist()/degree
print(phi[-1])
print((c.alpha + c.beta)/degree)
print((c.beta - c.alpha)/degree)

plt.subplot(3,1,1)
plt.plot(t, psi)
plt.ylabel(r'$\psi$  / deg')

plt.subplot(3,1,2)
plt.plot(t, gamma)
plt.ylabel(r'$\gamma$  / deg')

plt.subplot(3,1,3)
plt.plot(t, phi)
plt.ylabel(r'$\phi$  / deg')

plt.xlabel(r'$\theta$  / deg')
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
