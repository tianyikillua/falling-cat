# reference: T. R. Kane and M. P. Scher
#  "A Dynamical Explanation of the Falling Cat Phenomenon"
#    International Journal of Solids and Structures 5 (1969) 663

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp,cumtrapz
from mpl_toolkits.mplot3d import Axes3D

sqrt2 = np.sqrt(2)

class FallingCat:
    def __init__(self, JI, alpha, N=128, init=None):
        """
        inputs are:
        JI = (moment of intertia along cylinder axis)
            /(moment of intertia perpendicular to cylinder axis)
        alpha = angle between K and A1 in fig2 of Kane & Scher
        N = number of mesh points in theta [0,2pi]
        init = initial guess as an FallingCat object
        if init is not None, N is not used
        outputs as properties of FallingCat object are:
        beta = angle between K and B1 in fig2 of Kane & Scher
        theta = angle between A3 and B3
                as independent variable (shape(N,))
        psi = angle between spine plane and vertical plane
              as dependent variable (shape(N,))
        all angles are in radians
        """
        ca,sa = np.cos(alpha),np.sin(alpha)

        def fun(theta, psi, beta):
            """ eq(5) of Kane and Scher """
            cb,sb = np.cos(beta),np.sin(beta)
            ct = np.cos(theta)
            S = sqrt2*(ca*sb + sa*cb*ct)*sb
            T = ca*cb - sa*sb*ct
            T1,T2 = 1+T,1-T
            f = JI*S/T2/(T2 + JI*T1)/np.sqrt(T1)
            return [f]

        def bc(ya,yb,p):# boundary condition
            return np.r_[ya, yb-np.pi]

        if init is None:
            theta = np.linspace(0, 2*np.pi, N)
            psi = [np.linspace(0, np.pi, N)]
            beta = [(alpha + np.pi)/3]
        elif isinstance(init, FallingCat):
            theta = init.theta
            psi = [init.psi]
            beta = [init.beta]
        else:
            raise RuntimeError('bat init')

        s = solve_bvp(fun,bc,theta,psi,beta)

        self.alpha = alpha
        self.theta = s.x
        self.psi = s.y[0]
        self.beta = s.p[0]

    def bend(self, theta=None):
        """ gamma = angle between A1 and B1
                    in fig2 of Kane & Scher """
        if theta is None: theta = self.theta
        ca,sa = np.cos(self.alpha), np.sin(self.alpha)
        cb,sb = np.cos(self.beta), np.sin(self.beta)
        return np.arccos(ca*cb - sa*sb*np.cos(theta))

    def twist(self, theta=None):
        """ phi = angle of rotation along cylinder axis """
        a,b,t = self.alpha, self.beta, self.theta
        ca,sa = np.cos(a), np.sin(a)
        cb,sb = np.cos(b), np.sin(b)
        T11 = ca*cb - sa*sb*np.sin(t)
        T12 =-ca*sb - sa*cb*np.cos(t)
        p = cumtrapz(T12/(1 - T11**2), t, initial=0)*sb
        if theta is None: return p
        return np.interp(theta, t, p)

    def lean(self, theta=None):
        """ psi = angle between spine plane and vertical plane """
        if theta is None: return self.psi
        return np.interp(theta, self.theta, self.psi)

    def plot(self, theta, l=5, N=32,
             color=None, leg_color=None, **kw):
        """
        theta = angle between A3 and B3 in fig2 of Kane & Scher
        l = length of each cylinder
        N = number of mesh points on cylinder
        color = cylinder color (passed to Axes3D.plot_surface)
        leg_color = passed to Axes3D.plot
        kw = keyward arguments passed to Axes3D.plot
        assume Axes3D object has been already created
        """
        gam = self.bend(theta)
        phi = self.twist(theta)
        psi = self.lean(theta)
        cph,sph = np.cos(phi), np.sin(phi)
        cgm,sgm = np.cos(gam/2), np.sin(gam/2)
        cps,sps = np.cos(psi), np.sin(psi)

        ax = plt.gca()
        # plot cylinders
        u,v = np.mgrid[0:l:N*1j, 0:2*np.pi:N*1j]
        # front cylinder
        y = -u
        z,x = np.cos(v), np.sin(v)
        y,z = y*cgm + z*sgm,-y*sgm + z*cgm
        z,x = z*cps - x*sps, z*sps + x*cps
        ax.plot_surface(x,y,z,color=color)
        # rear cylinder
        y = u
        z,x = np.cos(v), np.sin(v)
        y,z = y*cgm - z*sgm, y*sgm + z*cgm
        z,x = z*cps - x*sps, z*sps + x*cps
        ax.plot_surface(x,y,z,color=color)
        # plot legs
        u = np.asarray([1,1])*l/2
        v = np.asarray([1,-1])*np.pi/12
        # front legs
        y = -u
        z,x = np.cos(v), np.sin(v)
        z,x = z*cph - x*sph, z*sph + x*cph
        y,z = y*cgm + z*sgm,-y*sgm + z*cgm
        z,x = z*cps - x*sps, z*sps + x*cps
        ax.plot(x,y,z,leg_color,**kw)
        # rear legs
        y = u
        z,x = np.cos(v), np.sin(v)
        z,x = z*cph - x*sph, z*sph + x*cph
        y,z = y*cgm - z*sgm, y*sgm + z*cgm
        z,x = z*cps - x*sps, z*sps + x*cps
        ax.plot(x,y,z,leg_color,**kw)

        return ax
