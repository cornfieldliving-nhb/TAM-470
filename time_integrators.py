import numpy as np
from pend_funcs import *

#Adams Bashforth 2.
def ab2pend(h, N, u0, ks, L0, m, g):
    t = np.linspace(0,N*h,N+1)
    urk2 = np.zeros((len(u0), N+1))
    urk2[:,0] = u0
    uab2 = np.zeros((len(u0), N+1))
    uab2[:,0] = u0
    alp = 1/2
    for i in range(N):
        ttem = t[i]
        utem = uab2[:,i]
        k1 = fpend(ttem, utem, ks, L0, m, g)
        k2 = fpend(ttem + alp*h, utem + alp*h*k1, ks, L0, m, g)
        urk2[:,i+1] = urk2[:,i] + (1 - 1/(2*alp))*h*k1 + 1/(2*alp)*h*k2
        if i == 0:
            uab2[:,i+1] = urk2[:,i+1]
        else:
            uab2[:,i+1] = uab2[:,i] + 1.5*h*fpend(ttem,utem,ks,L0,m,g) - 0.5*h*fpend(t[i-1],uab2[:,i-1],ks,L0,m,g)
    return uab2, t

#Leapfrog
def lfpend(h, N, u0, ks, L0, m, g):
    t = np.linspace(0,N*h,N+1)
    urk2 = np.zeros((len(u0), N+1))
    urk2[:,0] = u0
    ulf = np.zeros((len(u0), N+1))
    ulf[:,0] = u0
    alp = 1/2
    for i in range(N):
        ttem = t[i]
        utem = ulf[:,i]
        k1 = fpend(ttem, utem, ks, L0, m, g)
        k2 = fpend(ttem + alp*h, utem + alp*h*k1, ks, L0, m, g)
        urk2[:,i+1] = urk2[:,i] + (1 - 1/(2*alp))*h*k1 + 1/(2*alp)*h*k2
        if i == 0:
            ulf[:,i+1] = urk2[:,i+1]
        else:
            ulf[:,i+1] = ulf[:,i-1] + 2*h*fpend(ttem,utem,ks,L0,m,g)
    return ulf, t

#RK-4
def rk4pend(h, N, u0, ks, L0, m, g):
    t = np.linspace(0, N*h, N+1)
    urk4 = np.zeros((len(u0),N+1))
    urk4[:,0] = u0
    for i in range(N):
        ttem = t[i]
        utem = urk4[:,i]
        k1 = fpend(ttem, utem, ks, L0, m, g)
        k2 = fpend(ttem + h/2, utem + h/2*k1, ks, L0, m, g)
        k3 = fpend(ttem + h/2, utem + h/2*k2, ks, L0, m, g)
        k4 = fpend(ttem + h, utem + h*k3, ks, L0, m, g)
        urk4[:,i+1] = urk4[:,i] + h*(k1/6 + k2/3 + k3/3 + k4/6)
    return urk4,t

#TR-BDF2
def tr_bdf2pend(h, N, u0, ks, L0, m, g, tol=1e-10, max_iter=20):
    return #u, t