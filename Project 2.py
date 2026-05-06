import numpy as np
import matplotlib.pyplot as plt

#Fixed simulation parameters
Lx = 0.05
Ly = 0.03
w = 0.004
T_bar = 298
T_inf = 298
k = 48
rho = 7900
c = 470
sigma = 5.67e-8
epsilon = 0.8
beta = 2

def x_coord(vl,t):
    xlt = 0.1*Lx + vl*t
    return xlt

def laser_func(x, y, P, w, beta, Ly, xlt):
    r2 = (x-xlt)**2
    d = Ly - y
    f = beta/Ly * np.exp(-beta*(d/Ly))
    return (P/(w**2 *np.sqrt(2*np.pi))) * np.exp(-r2/(2*w**2)) * f

#Problem 2 
P_2 = 800
Nx_2 = 50
Ny_2 = 50
tend_2 = 150
x_2 = Lx/2
y_2 = Ly

x_2_field = np.linspace(0,Lx,Nx_2+1)
y_2_field = np.linspace(0,Ly,Ny_2+1)
dt = 0.01
dx = Lx/(Nx_2)
dy = Ly/(Ny_2)
Nt_2 = int(tend_2/dt)
t_2 = np.linspace(0,tend_2,Nt_2+1)
T = np.ones((Nx_2+1, Ny_2+1, Nt_2+1)) * 298 




conv_cri = 1e-6
max_ite = 50

for tind in range(Nt_2):
    for i in range(1,Nx_2):
        for j in range(1,Ny_2):
            D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
            D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
            T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_2_field[i], y_2_field[j], P_2, w, beta, Ly, x_2)
    
    
    #left-side BC
    for y_ind in range(Ny_2+1):
        T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
        for ite1 in range(max_ite):
            f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
            f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
            delta = f_nt_ls/f_d_nt_ls
            T_old = T[0,y_ind,tind+1]    
            T[0,y_ind,tind+1] -= delta
            if (abs(f_nt_ls)) < conv_cri:
                break
            # print("ite1:" + str(ite1))

    # #right-side BC
    for y_ind2 in range(Ny_2+1):
        T[Nx_2,y_ind2,tind+1] = T[Nx_2-1,y_ind2,tind+1]
        for ite2 in range(max_ite):
            f_nt_rs = (T[Nx_2,y_ind2,tind+1] - T[Nx_2-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_2,y_ind2,tind+1])**4 - T_inf**4)/k
            f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_2,y_ind2,tind+1])**3)
            delta = f_nt_rs/f_d_nt_rs
            T_old = T[Nx_2,y_ind2,tind+1]
            T[Nx_2,y_ind2,tind+1] -= delta
            if (abs(f_nt_rs)) < conv_cri:
                break
            # print("ite2:" + str(ite2))
        
    # #top-side BC
    for x_ind in range(Nx_2+1):
        T[x_ind,Ny_2,tind+1] = T[x_ind,Ny_2-1,tind+1]
        for ite3 in range(max_ite):
            f_nt_ts = (T[x_ind,Ny_2,tind+1] - T[x_ind,Ny_2-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_2,tind+1])**4 - T_inf**4)/k
            f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_2,tind+1])**3)
            delta = f_nt_ts/f_d_nt_ts
            T_old = T[x_ind,Ny_2,tind+1]
            T[x_ind,Ny_2,tind+1] -= delta
            if (abs(f_nt_ts)) < conv_cri:
                break
            # print("ite3:" + str(ite3) )
        
    # #bottom-side BC    
    T[:,0,tind+1] = T_bar

#Problem 2 plot
#P2-1
times = [30, 60, 90, 120, 150]
print("\nMaximum Temperature at Top Surface:")
print("Time (s) | Temperature (K)")
print("-" * 30)
for t in times:
    t_index = int(t/dt)
    max_temp = np.max(T[:, Ny_2, t_index])  
    print(f"{t:8.0f} | {max_temp:13.2f}")


#p2-2
t_2 = np.linspace(0, tend_2, Nt_2+1)


x_laser = int(Nx_2/2)           
y_top = Ny_2                 
y_20 = int(0.8*Ny_2)         
y_50 = int(0.5*Ny_2)         

# Get temperature histories at these points
T_top = T[x_laser, y_top, :]    # Temperature at top surface
T_20 = T[x_laser, y_20, :]      # Temperature at 20% depth
T_50 = T[x_laser, y_50, :]      # Temperature at 50% depth

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(t_2, T_top, 'r-', label='Top surface')
plt.plot(t_2, T_20, 'g-', label='20% depth')
plt.plot(t_2, T_50, 'b-', label='50% depth')

plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature curve at different depths (x = 0.5Lx)')
plt.legend()
plt.grid(True)
plt.show()

#p2-3
x_mid = int(Nx_2/2)  # x = 0.5Lx
T_profile = T[x_mid, :, -1]  
y_positions = y_2_field 
plt.figure(figsize=(10, 6))
plt.plot(y_positions, T_profile, 'b-', linewidth=2)

plt.xlabel('y position (m)')
plt.ylabel('Temperature (K)')
plt.title(f'Temperature vs Position at x = 0.5Lx, t = {tend_2}s')
plt.grid(True)
plt.show()

#p2-4
times = [10, 40, 150]
for t in times:
    t_index = int(t/dt)
    plt.figure(figsize = (10,6))
    im = plt.contourf(x_2_field,y_2_field,T[:,:,t_index].T,levels=20,cmap='rainbow')
    plt.colorbar(im, label='Temperature (K)')
    plt.xlabel('x(m)')
    plt.ylabel('y(m)')
    plt.title(f'temperature Distribution at t = {t}s')
    plt.show()


#Problem 3
#p3-1
P_2 = 800
vl = 0.01
x_start = 0.1*Lx
x_end = 0.9*Lx
Nx_3 = 50
Ny_3 = 50
tend_3 = (x_end - x_start)/vl


x_3_field = np.linspace(0,Lx,Nx_3+1)
y_3_field = np.linspace(0,Ly,Ny_3+1)
dt = 0.01
dx = Lx/(Nx_3)
dy = Ly/(Ny_3)
Nt_3 = int(tend_3/dt)
t_3 = np.linspace(0,tend_3,Nt_3+1)
T = np.ones((Nx_3+1, Ny_3+1, Nt_3+1)) * 298 

conv_cri = 1e-6
max_ite = 50

for tind in range(Nt_3):
    current_time = tind*dt
    x = x_coord(vl,current_time)
    for i in range(1,Nx_3):
        for j in range(1,Ny_3):
            D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
            D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
            T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], P_2, w, beta, Ly, x)
    
    
    #left-side BC
    for y_ind in range(Ny_3+1):
        T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
        for ite1 in range(max_ite):
            f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
            f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
            delta = f_nt_ls/f_d_nt_ls
            T_old = T[0,y_ind,tind+1]    
            T[0,y_ind,tind+1] -= delta
            if (abs(f_nt_ls)) < conv_cri:
                break
            # print("ite1:" + str(ite1))

    # #right-side BC
    for y_ind2 in range(Ny_3+1):
        T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
        for ite2 in range(max_ite):
            f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
            f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
            delta = f_nt_rs/f_d_nt_rs
            T_old = T[Nx_3,y_ind2,tind+1]
            T[Nx_3,y_ind2,tind+1] -= delta
            if (abs(f_nt_rs)) < conv_cri:
                break
            # print("ite2:" + str(ite2))
        
    # #top-side BC
    for x_ind in range(Nx_3+1):
        T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
        for ite3 in range(max_ite):
            f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
            f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
            delta = f_nt_ts/f_d_nt_ts
            T_old = T[x_ind,Ny_3,tind+1]
            T[x_ind,Ny_3,tind+1] -= delta
            if (abs(f_nt_ts)) < conv_cri:
                break
            # print("ite3:" + str(ite3) )
        
    # #bottom-side BC    
    T[:,0,tind+1] = T_bar

t_3 = np.linspace(0, tend_3, Nt_3+1)


x_laser = int(Nx_3/2)           
y_top = Ny_2                 
y_20 = int(0.8*Ny_2)         
y_50 = int(0.5*Ny_2)         

# Get temperature histories at these points
T_top = T[x_laser, y_top, :]    # Temperature at top surface
T_20 = T[x_laser, y_20, :]      # Temperature at 20% depth
T_50 = T[x_laser, y_50, :]      # Temperature at 50% depth

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(t_3, T_top, 'r-', label='Top surface')
plt.plot(t_3, T_20, 'g-', label='20% depth')
plt.plot(t_3, T_50, 'b-', label='50% depth')

plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature curve at different depths (x = 0.5Lx)')
plt.legend()
plt.grid(True)
plt.show()



#p3-3-1
#t = 0.2*tend
for tind in range(Nt_3):
    current_time = tind*dt
    x_3 = x_coord(vl,current_time)
    for i in range(1,Nx_3):
        for j in range(1,Ny_3):
            D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
            D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
            T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], P_2, w, beta, Ly, x_3)
    
    
    #left-side BC
    for y_ind in range(Ny_3+1):
        T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
        for ite1 in range(max_ite):
            f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
            f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
            delta = f_nt_ls/f_d_nt_ls
            T_old = T[0,y_ind,tind+1]    
            T[0,y_ind,tind+1] -= delta
            if (abs(f_nt_ls)) < conv_cri:
                break
            # print("ite1:" + str(ite1))

    # #right-side BC
    for y_ind2 in range(Ny_3+1):
        T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
        for ite2 in range(max_ite):
            f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
            f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
            delta = f_nt_rs/f_d_nt_rs
            T_old = T[Nx_3,y_ind2,tind+1]
            T[Nx_3,y_ind2,tind+1] -= delta
            if (abs(f_nt_rs)) < conv_cri:
                break
            # print("ite2:" + str(ite2))
        
    # #top-side BC
    for x_ind in range(Nx_3+1):
        T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
        for ite3 in range(max_ite):
            f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
            f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
            delta = f_nt_ts/f_d_nt_ts
            T_old = T[x_ind,Ny_3,tind+1]
            T[x_ind,Ny_3,tind+1] -= delta
            if (abs(f_nt_ts)) < conv_cri:
                break
            # print("ite3:" + str(ite3) )
        
    # #bottom-side BC    
    T[:,0,tind+1] = T_bar

plt.figure(figsize = (10,6))
t_index = int(0.2*tend_3/dt)
im = plt.contourf(x_3_field, y_3_field, T[:,:,t_index].T, levels=20, cmap='rainbow')  
plt.colorbar(im, label='Temperature (K)')
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.title(f'temperature Distribution at t = {0.2*tend_3}s')
plt.show()

#p3-3-2
#t = 0.5*tend
for tind in range(Nt_3):
    current_time = tind*dt
    x_3 = x_coord(vl,current_time)
    for i in range(1,Nx_3):
        for j in range(1,Ny_3):
            D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
            D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
            T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], P_2, w, beta, Ly, x_3)
    
    
    #left-side BC
    for y_ind in range(Ny_3+1):
        T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
        for ite1 in range(max_ite):
            f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
            f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
            delta = f_nt_ls/f_d_nt_ls
            T_old = T[0,y_ind,tind+1]    
            T[0,y_ind,tind+1] -= delta
            if (abs(f_nt_ls)) < conv_cri:
                break
            # print("ite1:" + str(ite1))

    # #right-side BC
    for y_ind2 in range(Ny_3+1):
        T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
        for ite2 in range(max_ite):
            f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
            f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
            delta = f_nt_rs/f_d_nt_rs
            T_old = T[Nx_3,y_ind2,tind+1]
            T[Nx_3,y_ind2,tind+1] -= delta
            if (abs(f_nt_rs)) < conv_cri:
                break
            # print("ite2:" + str(ite2))
        
    # #top-side BC
    for x_ind in range(Nx_3+1):
        T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
        for ite3 in range(max_ite):
            f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
            f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
            delta = f_nt_ts/f_d_nt_ts
            T_old = T[x_ind,Ny_3,tind+1]
            T[x_ind,Ny_3,tind+1] -= delta
            if (abs(f_nt_ts)) < conv_cri:
                break
            # print("ite3:" + str(ite3) )
        
    # #bottom-side BC    
    T[:,0,tind+1] = T_bar

plt.figure(figsize = (10,6))
t_index = int(0.5*tend_3/dt)
im = plt.contourf(x_3_field,y_3_field,T[:,:,t_index].T,levels=20,cmap='rainbow')
plt.colorbar(im, label='Temperature (K)')
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.title(f'temperature Distribution at t = {0.5*tend_3}s')
plt.show()


#p3-3-3
#t = tend
for tind in range(Nt_3):
    current_time = tind*dt
    x_3 = x_coord(vl,current_time)
    for i in range(1,Nx_3):
        for j in range(1,Ny_3):
            D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
            D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
            T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], P_2, w, beta, Ly, x_3)
    
    
    #left-side BC
    for y_ind in range(Ny_3+1):
        T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
        for ite1 in range(max_ite):
            f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
            f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
            delta = f_nt_ls/f_d_nt_ls
            T_old = T[0,y_ind,tind+1]    
            T[0,y_ind,tind+1] -= delta
            if (abs(f_nt_ls)) < conv_cri:
                break
            # print("ite1:" + str(ite1))

    # #right-side BC
    for y_ind2 in range(Ny_3+1):
        T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
        for ite2 in range(max_ite):
            f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
            f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
            delta = f_nt_rs/f_d_nt_rs
            T_old = T[Nx_3,y_ind2,tind+1]
            T[Nx_3,y_ind2,tind+1] -= delta
            if (abs(f_nt_rs)) < conv_cri:
                break
            # print("ite2:" + str(ite2))
        
    # #top-side BC
    for x_ind in range(Nx_3+1):
        T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
        for ite3 in range(max_ite):
            f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
            f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
            delta = f_nt_ts/f_d_nt_ts
            T_old = T[x_ind,Ny_3,tind+1]
            T[x_ind,Ny_3,tind+1] -= delta
            if (abs(f_nt_ts)) < conv_cri:
                break
            # print("ite3:" + str(ite3) )
        
    # #bottom-side BC    
    T[:,0,tind+1] = T_bar

plt.figure(figsize = (10,6))
t_index = int(tend_3/dt)
im = plt.contourf(x_3_field,y_3_field,T[:,:,t_index].T,levels=20,cmap='rainbow')
plt.colorbar(im, label='Temperature (K)')
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.title(f'temperature Distribution at t = {tend_3}s')
plt.show()


#Problem 4
#p4-1
P = [200, 400, 600, 800, 1000, 1200]
vl_mat = [0.01,0.02,0.04,0.08]
#vl = 0.01
x_start = 0.1*Lx
x_end = 0.9*Lx
Nx_3 = 50
Ny_3 = 50
tend_3 = (x_end - x_start)/vl


x_3_field = np.linspace(0,Lx,Nx_3+1)
y_3_field = np.linspace(0,Ly,Ny_3+1)
dt = 0.01
dx = Lx/(Nx_3)
dy = Ly/(Ny_3)
Nt_3 = int(tend_3/dt)
t_3 = np.linspace(0,tend_3,Nt_3+1)
T = np.ones((Nx_3+1, Ny_3+1, Nt_3+1)) * 298 

conv_cri = 1e-6
max_ite = 50
T_val_v1 = []
T_val_v2 = []
T_val_v3 = []
T_val_v4 = []

#V1
for p_val in P:
    for tind in range(Nt_3):
        current_time = tind*dt
        x = x_coord(vl_mat[0],current_time)
        for i in range(1,Nx_3):
            for j in range(1,Ny_3):
                D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
                D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
                T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], p_val, w, beta, Ly, x)
        
        
        #left-side BC
        for y_ind in range(Ny_3+1):
            T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
            for ite1 in range(max_ite):
                f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
                f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
                delta = f_nt_ls/f_d_nt_ls
                T_old = T[0,y_ind,tind+1]    
                T[0,y_ind,tind+1] -= delta
                if (abs(f_nt_ls)) < conv_cri:
                    break
                # print("ite1:" + str(ite1))

        # #right-side BC
        for y_ind2 in range(Ny_3+1):
            T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
            for ite2 in range(max_ite):
                f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
                f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
                delta = f_nt_rs/f_d_nt_rs
                T_old = T[Nx_3,y_ind2,tind+1]
                T[Nx_3,y_ind2,tind+1] -= delta
                if (abs(f_nt_rs)) < conv_cri:
                    break
                # print("ite2:" + str(ite2))
            
        # #top-side BC
        for x_ind in range(Nx_3+1):
            T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
            for ite3 in range(max_ite):
                f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
                f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
                delta = f_nt_ts/f_d_nt_ts
                T_old = T[x_ind,Ny_3,tind+1]
                T[x_ind,Ny_3,tind+1] -= delta
                if (abs(f_nt_ts)) < conv_cri:
                    break
                # print("ite3:" + str(ite3) )
            
        # #bottom-side BC    
        T[:,0,tind+1] = T_bar
    T_val_v1.append(np.max(T))

#V2
for p_val in P:
    for tind in range(Nt_3):
        current_time = tind*dt
        x = x_coord(vl_mat[1],current_time)
        for i in range(1,Nx_3):
            for j in range(1,Ny_3):
                D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
                D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
                T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], p_val, w, beta, Ly, x)
        
        
        #left-side BC
        for y_ind in range(Ny_3+1):
            T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
            for ite1 in range(max_ite):
                f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
                f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
                delta = f_nt_ls/f_d_nt_ls
                T_old = T[0,y_ind,tind+1]    
                T[0,y_ind,tind+1] -= delta
                if (abs(f_nt_ls)) < conv_cri:
                    break
                # print("ite1:" + str(ite1))

        # #right-side BC
        for y_ind2 in range(Ny_3+1):
            T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
            for ite2 in range(max_ite):
                f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
                f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
                delta = f_nt_rs/f_d_nt_rs
                T_old = T[Nx_3,y_ind2,tind+1]
                T[Nx_3,y_ind2,tind+1] -= delta
                if (abs(f_nt_rs)) < conv_cri:
                    break
                # print("ite2:" + str(ite2))
            
        # #top-side BC
        for x_ind in range(Nx_3+1):
            T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
            for ite3 in range(max_ite):
                f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
                f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
                delta = f_nt_ts/f_d_nt_ts
                T_old = T[x_ind,Ny_3,tind+1]
                T[x_ind,Ny_3,tind+1] -= delta
                if (abs(f_nt_ts)) < conv_cri:
                    break
                # print("ite3:" + str(ite3) )
            
        # #bottom-side BC    
        T[:,0,tind+1] = T_bar
    T_val_v2.append(np.max(T))

#V3
for p_val in P:
    for tind in range(Nt_3):
        current_time = tind*dt
        x = x_coord(vl_mat[2],current_time)
        for i in range(1,Nx_3):
            for j in range(1,Ny_3):
                D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
                D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
                T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], p_val, w, beta, Ly, x)
        
        
        #left-side BC
        for y_ind in range(Ny_3+1):
            T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
            for ite1 in range(max_ite):
                f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
                f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
                delta = f_nt_ls/f_d_nt_ls
                T_old = T[0,y_ind,tind+1]    
                T[0,y_ind,tind+1] -= delta
                if (abs(f_nt_ls)) < conv_cri:
                    break
                # print("ite1:" + str(ite1))

        # #right-side BC
        for y_ind2 in range(Ny_3+1):
            T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
            for ite2 in range(max_ite):
                f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
                f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
                delta = f_nt_rs/f_d_nt_rs
                T_old = T[Nx_3,y_ind2,tind+1]
                T[Nx_3,y_ind2,tind+1] -= delta
                if (abs(f_nt_rs)) < conv_cri:
                    break
                # print("ite2:" + str(ite2))
            
        # #top-side BC
        for x_ind in range(Nx_3+1):
            T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
            for ite3 in range(max_ite):
                f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
                f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
                delta = f_nt_ts/f_d_nt_ts
                T_old = T[x_ind,Ny_3,tind+1]
                T[x_ind,Ny_3,tind+1] -= delta
                if (abs(f_nt_ts)) < conv_cri:
                    break
                # print("ite3:" + str(ite3) )
            
        # #bottom-side BC    
        T[:,0,tind+1] = T_bar
    T_val_v3.append(np.max(T))

#V4
for p_val in P:
    for tind in range(Nt_3):
        current_time = tind*dt
        x = x_coord(vl_mat[3],current_time)
        for i in range(1,Nx_3):
            for j in range(1,Ny_3):
                D2x = (T[i+1,j,tind] - 2*T[i,j,tind] + T[i-1,j,tind])/(dx)**2
                D2y = (T[i,j+1,tind] - 2*T[i,j,tind] + T[i,j-1,tind])/(dy)**2
                T[i,j,tind+1] = T[i,j,tind] + k*dt/(rho*c)*(D2x + D2y) + dt/(rho*c) * laser_func(x_3_field[i], y_3_field[j], p_val, w, beta, Ly, x)
        
        
        #left-side BC
        for y_ind in range(Ny_3+1):
            T[0,y_ind,tind+1] = T[1,y_ind,tind+1]
            for ite1 in range(max_ite):
                f_nt_ls = (T[0,y_ind,tind+1] - T[1,y_ind,tind+1])*(1/dx) + sigma*epsilon*((T[0,y_ind,tind+1])**4 - T_inf**4)/k
                f_d_nt_ls = 1/dx + sigma*epsilon/k * (4*(T[0,y_ind,tind+1])**3)
                delta = f_nt_ls/f_d_nt_ls
                T_old = T[0,y_ind,tind+1]    
                T[0,y_ind,tind+1] -= delta
                if (abs(f_nt_ls)) < conv_cri:
                    break
                # print("ite1:" + str(ite1))

        # #right-side BC
        for y_ind2 in range(Ny_3+1):
            T[Nx_3,y_ind2,tind+1] = T[Nx_3-1,y_ind2,tind+1]
            for ite2 in range(max_ite):
                f_nt_rs = (T[Nx_3,y_ind2,tind+1] - T[Nx_3-1,y_ind2,tind+1])*(1/dx) + sigma*epsilon*((T[Nx_3,y_ind2,tind+1])**4 - T_inf**4)/k
                f_d_nt_rs = 1/dx + sigma*epsilon/k * (4*(T[Nx_3,y_ind2,tind+1])**3)
                delta = f_nt_rs/f_d_nt_rs
                T_old = T[Nx_3,y_ind2,tind+1]
                T[Nx_3,y_ind2,tind+1] -= delta
                if (abs(f_nt_rs)) < conv_cri:
                    break
                # print("ite2:" + str(ite2))
            
        # #top-side BC
        for x_ind in range(Nx_3+1):
            T[x_ind,Ny_3,tind+1] = T[x_ind,Ny_3-1,tind+1]
            for ite3 in range(max_ite):
                f_nt_ts = (T[x_ind,Ny_3,tind+1] - T[x_ind,Ny_3-1,tind+1])*(1/dy) + sigma*epsilon*((T[x_ind,Ny_3,tind+1])**4 - T_inf**4)/k
                f_d_nt_ts = 1/dy + sigma*epsilon/k * (4*(T[x_ind,Ny_3,tind+1])**3)
                delta = f_nt_ts/f_d_nt_ts
                T_old = T[x_ind,Ny_3,tind+1]
                T[x_ind,Ny_3,tind+1] -= delta
                if (abs(f_nt_ts)) < conv_cri:
                    break
                # print("ite3:" + str(ite3) )
            
        # #bottom-side BC    
        T[:,0,tind+1] = T_bar
    T_val_v4.append(np.max(T))

plt.figure(figsize = (10,6))
plt.plot(P, T_val_v1, label='v1 condition')
plt.plot(P, T_val_v2, label='v2 condition')
plt.plot(P, T_val_v3, label='v3 condition')
plt.plot(P, T_val_v4, label='v4 condition')
plt.xlabel('power (W)')
plt.ylabel('Temp (K)')
plt.grid()
plt.legend()
plt.show()