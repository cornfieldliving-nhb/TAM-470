import numpy as np

#The function which generates the f vector.
def fpend(t,u,ks,L0,m,g):
    f = np.zeros(len(u))
    f[0:4] = u[4:len(u)]
    L1 = np.sqrt(u[0]**2 + u[1]**2)
    L2 = np.sqrt((u[2]-u[0])**2 + (u[3]-u[1])**2)
    f[-4] = (-ks[0]*(L1 - L0[0])* u[0]/L1 + ks[1]*(L2 - L0[1])*(u[2]-u[0])/L2)/m[0]
    f[-3] = (ks[0]*(L1 - L0[0])* -u[1]/L1 - ks[1]*(L2 - L0[1])*(-(u[3]-u[1]))/L2 - m[0]*g)/m[0]
    f[-2] = (-ks[1] * (L2 - L0[1]) * (u[2]-u[0])/L2)/m[1]
    f[-1] = (ks[1] * (L2 - L0[1]) * (-(u[3] - u[1]))/L2 - m[1]*g)/m[1]
    return f




#Gradient of fpend.
def dfpend(t,u,ks,L0,m,g):
    return #df

def pendulumEnergy(U,ks,L0,m,g):
    E = np.zeros(U.shape[1])
    for i in range(U.shape[1]):
        utem = U[:,i]
        L1 = np.sqrt(utem[0]**2 + utem[1]**2)
        L2 = np.sqrt((utem[2]-utem[0])**2 + (utem[3]-utem[1])**2)
        KE = 1/2*m[0]*(utem[4]**2 + utem[5]**2) + 1/2*m[1]*(utem[6]**2 + utem[7]**2)
        PE = 1/2*ks[0]*(L1-L0[0])**2 + 1/2*ks[1]*(L2-L0[1])**2 + m[0]*g*utem[1] + m[1]*g*utem[3]
        E[i] += KE + PE
    return E