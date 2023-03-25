import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import degree
from FallingCat import FallingCat

JI = [0.2, 0.3]
alpha = np.linspace(36,82,50)*degree

plt.figure(figsize=(5,3.75))

for j in JI:
    beta,c = [],None
    for a in alpha:
        c = FallingCat(j,a,init=c)
        beta.append(c.beta)

    ba = (beta - alpha)/degree
    ab = (alpha + beta)/degree
    plt.plot(ba, ab, label='I = %g'%j)

plt.legend()
plt.xlabel(r'$\beta - \alpha$  / deg')
plt.ylabel(r'$\alpha + \beta$  / deg')
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
