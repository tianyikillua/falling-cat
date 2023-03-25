import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import degree
from FallingCat import FallingCat

JI = 0.25
alpha = 30*degree
theta = np.linspace(0, 2*np.pi, 5)

c = FallingCat(JI, alpha)

print(c.beta/degree)
print(c.bend(theta))
print(c.twist(theta))
print(c.lean(theta))

plt.figure(figsize=(2,6))

for i,th in enumerate(theta):
    plt.subplot(5,1,i+1, projection='3d')
    ax = c.plot(th, color='y', leg_color='.b', ms=10)
    ax.view_init(15, -15)
    ax.set_axis_off()
    ax.autoscale(False)
#    ax.set_xlim(-4,4)
#    ax.set_ylim(-4,4)
#    ax.set_zlim(-4,4)
#    ax.axis('equal')

plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
