import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
from time_integrators import *
from pend_funcs import *
from utilities import *

def main():
    '''
    ########To use scipy_ivp, use the following as a sample code:#############
    ks = np.array([1e6,1e6]) #Stifness values of both springs.
    L0 = np.array([1,1])
    m = np.array([1,1]) #Mass of the two particles.

    #Simulation time details
    Tf = 600 #Final time.
    h = 0.01 #Time step size.
    N = int(Tf/h) #Number of steps based on above parameters.
    g = 9.81
    t_span = [0,Tf]

    #IC
    u0 = np.array([1,0,1,1,0,0,0,0])

    #Scipy solution
    sol = spi.solve_ivp(fpend, t_span, u0, method='BDF',args=(ks,L0,m,g))
    t = sol.t; y = sol.y
    #####################################################################################

    #########################Instructions for animation##############################################
    #See below for instructions to generate animations. May not work in the workspace.
    #Need to install ffmpeg using conda or pip. - (conda install -c conda-forge ffmpeg) or (pip install ffmpeg OR pip3 install ffmpeg)

    #There are three ways you can generate an animation:-

    #1. When you want to use every time step to generate the animation. This will be extremely slow when there are a lot of time steps.
    #Useful when you want to see the animation at every time step and make sure that t.shape[0] is a smaller number (say below 5000)
    #pend_plot_gen_animation(t,y,animation_length=10, filename="pend.mp4")

    #2. When you want to generate an animation of fixed length and fixed fps. This will be quick depending on the parametrs but may not be very #useful to interpret.
    #pend_plot_gen_animation_fast(t, y, animation_length=10, fps=60, filename="pend_fast.mp4")

    #3. When you want an animation of only the first few time steps with a specified fps. This can be useful to see if the system is behaving well #initially and can be pretty fast depending on the parameter n. This is preferred when t.shape[0] >>> n.
    #Define a parameter n below to choose how many initial time steps you want to see.
    #pend_plot_gen_animation_first_n_frames(t, y, n=1200, fps=60, filename="pend_first_n.mp4")
    ###########################################################################
    '''

    ks = np.array([500,500]) #Stifness values of both springs.
    L0 = np.array([1,1])
    m = np.array([1,1]) #Mass of the two particles.

    #Simulation time details
    Tf = 600 #Final time.
    h = 0.01 #Time step size.
    N = int(Tf/h) #Number of steps based on above parameters.
    g = 9.81
    t_span = [0,Tf]

    #IC
    u0 = np.array([1,0,1,1,0,0,0,0])

    #print(fpend(0,u0,ks,L0,m,g))
    # sol = spi.solve_ivp(fpend, t_span, u0, method='Radau',args=(ks,L0,m,g))
    # t = sol.t; y = sol.y
    # Esci = pendulumEnergy(y,ks,L0,m,g)
    #pend_plot_gen_animation_first_n_frames(t, y, n=t.shape[0], fps=60, filename="scipend_first_n.mp4")
    

    #h=0.02
    h1=0.02
    N1 = int(Tf/h1)
    Urk4_1, trk4_1 = rk4pend(h1, N1, u0, ks, L0, m, g)
    Erk4_1 = pendulumEnergy(Urk4_1,ks,L0,m,g)

    Ulf_1,tlf_1 = lfpend(h1,N1,u0,ks,L0,m,g)
    Elf_1 = pendulumEnergy(Ulf_1,ks,L0,m,g)

    Uab2_1,tab2_1 = ab2pend(h1,N1,u0,ks,L0,m,g)
    Eab2_1 = pendulumEnergy(Uab2_1,ks,L0,m,g)

    inte_size = 10
    i = 0
    inte_avg_lf1 = []
    for i in range(len(Elf_1)):
        tem_inte = Elf_1[i:i+inte_size]
        tem_inte = np.average(Elf_1)
        inte_avg_lf1.append(tem_inte)
        i += inte_size


    # plt.figure("Energy VS time with h=0.02")
    # plt.title("Energy VS time with h=0.02")
    # plt.subplot(1,3,1)
    # plt.plot(trk4_1,Erk4_1,'--',color = 'blue')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('rk4')
    # plt.subplot(1,3,2)
    # plt.plot(tlf_1,Elf_1,'--',color = 'red')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('leap frog')
    # plt.subplot(1,3,3)
    # plt.plot(tab2_1,Eab2_1,'--',color = 'green')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('ab2')
    # plt.show()

    # plt.figure("Energy VS time with h=0.01")
    # plt.title("Energy VS time with h=0.01")
    # plt.subplot(1,3,1)
    # plt.plot(trk4_2,Erk4_2,'--',color = 'blue')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('rk4')
    # plt.subplot(1,3,2)
    # plt.plot(tlf_2,Elf_2,'--',color = 'red')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('leap frog')
    # plt.subplot(1,3,3)
    # plt.plot(tab2_2,Eab2_2,'--',color = 'green')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('ab2')
    # plt.show()

    # plt.figure("Energy VS time with h=0.005")
    # plt.title("Energy VS time with h=0.005")
    # plt.subplot(1,3,1)
    # plt.plot(trk4_3,Erk4_3,'--',color = 'blue')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('rk4')
    # plt.subplot(1,3,2)
    # plt.plot(tlf_3,Elf_3,'--',color = 'red')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('leap frog')
    # plt.subplot(1,3,3)
    # plt.plot(tab2_3,Eab2_3,'--',color = 'green')
    # plt.xlabel("t")
    # plt.ylabel("Energy")
    # plt.title('ab2')
    # plt.show()

    # pend_plot_gen_animation_first_n_frames(trk4_1, Urk4_1, n=t.shape[0], fps=60, filename="scipend_first_n.mp4")
    # pend_plot_gen_animation_first_n_frames(tlf_1, Ulf_1, n=t.shape[0], fps=60, filename="scipend_first_n.mp4")
    # pend_plot_gen_animation_first_n_frames(tab2_1, Uab2_1, n=t.shape[0], fps=60, filename="scipend_first_n.mp4")

    #h=0.01
    h2 = 0.01
    N2 = int(Tf/h2)
    Urk4_2, trk4_2 = rk4pend(h2, N2, u0, ks, L0, m, g)
    Erk4_2 = pendulumEnergy(Urk4_2,ks,L0,m,g)

    Ulf_2,tlf_2 = lfpend(h2,N2,u0,ks,L0,m,g)
    Elf_2 = pendulumEnergy(Ulf_2,ks,L0,m,g)

    Uab2_2,tab2_2 = ab2pend(h2,N2,u0,ks,L0,m,g)
    Eab2_2 = pendulumEnergy(Uab2_2,ks,L0,m,g)

    inte_size = 10
    i = 0
    inte_avg_lf2 = []
    for i in range(len(Elf_2)):
        tem_inte = Elf_2[i:i+inte_size]
        tem_inte = np.average(Elf_2)
        inte_avg_lf2.append(tem_inte)
        i += inte_size


    
    #h=0.005
    h3 = 0.005
    N3 = int(Tf/h3)
    Urk4_3, trk4_3 = rk4pend(h3, N3, u0, ks, L0, m, g)
    Erk4_3 = pendulumEnergy(Urk4_3,ks,L0,m,g)

    Ulf_3,tlf_3 = lfpend(h3,N3,u0,ks,L0,m,g)
    Elf_3 = pendulumEnergy(Ulf_3,ks,L0,m,g)

    Uab2_3,tab2_3 = ab2pend(h3,N3,u0,ks,L0,m,g)
    Eab2_3 = pendulumEnergy(Uab2_3,ks,L0,m,g)

    inte_size = 10
    i = 0
    inte_avg_lf3 = []
    for i in range(len(Elf_3)):
        tem_inte = Elf_3[i:i+inte_size]
        tem_inte = np.average(Elf_3)
        inte_avg_lf3.append(tem_inte)
        i += inte_size

    plt.figure("Average Total Energy vs time")
    plt.plot(tlf_1, inte_avg_lf1,'--',label = "lf h=0.02")
    plt.plot(tlf_2, inte_avg_lf2,'--',label = "lf h=0.01")
    plt.plot(tlf_3, inte_avg_lf3,'--',label = "lf h=0.005")
    plt.plot(trk4_1, Erk4_1,'--',label = "rk4 h=0.02")
    plt.plot(trk4_2, Erk4_2,'--',label = "rk4 h=0.01")
    plt.plot(trk4_3, Erk4_3,'--',label = "rk4 h=0.005")
    plt.legend()
    plt.xlabel("t")
    plt.ylabel("Avg E")
    plt.show()

    # print("Lf average energy at h=0.02 is:" + str(np.average(inte_avg_lf2)))
    # print("LF average energy at h=0.01 is:" + str(np.average(inte_avg_lf2)))
    # print("LF average energy at h=0.005 is:" + str(np.average(inte_avg_lf3)))

    
    
    
    #Question 3
    ks_2 = np.array([1e6,1e6]) #Stifness values of both springs.
    L0 = np.array([1,1])
    m = np.array([1,1]) #Mass of the two particles.

    #Simulation time details
    Tf = 100 #Final time.
    h = 0.01 #Time step size.
    N = int(Tf/h) #Number of steps based on above parameters.
    g = 9.81
    t_span = [0,Tf]

    #IC
    u0_2 = np.array([1,0,1,1,0,0,0,0])

    #print(fpend(0,u0,ks,L0,m,g))
    sol = spi.solve_ivp(fpend, t_span, u0_2, method='Radau',args=(ks_2,L0,m,g))
    t = sol.t; y = sol.y
    Esci = pendulumEnergy(y,ks_2,L0,m,g)

    plt.figure("Radau Energy vs time")
    plt.plot(t,Esci,'--')
    plt.xlabel("t")
    plt.ylabel("Total energy")
    plt.show()

    #pend_plot_gen_animation_first_n_frames(t, y, n=t.shape[0], fps=60, filename="scipend_radau_n.mp4")

    sol_2 = spi.solve_ivp(fpend, t_span, u0_2, method='RK45',args=(ks_2,L0,m,g))
    t_2 = sol_2.t; y_2 = sol_2.y
    Esci_2 = pendulumEnergy(y_2,ks_2,L0,m,g)

    plt.figure("RK45 Energy vs time")
    plt.plot(t_2,Esci_2,'--')
    plt.xlabel("t")
    plt.ylabel("Total energy")
    plt.show()

    #pend_plot_gen_animation_first_n_frames(t_2, y_2, n=t.shape[0], fps=60, filename="scipend_rk45_n.mp4")

    
    


    
    



#No need to change anything below this line.
if __name__ == "__main__":
    main()