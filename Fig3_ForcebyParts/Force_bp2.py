# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 18:24:00 2016

@author: Enrique Alejandro
"""


from flatpunch import z_eq, force_by_parts, figannotate
import numpy as np
import matplotlib.pyplot as plt



a = 0.2e-9  #Minimum intermolecular distance
params_simul = np.loadtxt('Summary.txt', skiprows=1)
Zr = params_simul[7]*1.0e-9
R = params_simul[3]*1.0e-9  #radius of the tip
nu = params_simul[4]   #Poison ratio
H = params_simul[5]*1.0e-19  #Hammaker constant
alfa = 4.0*R/(1.0-nu) #Cell constant converting stress/strain to force/displacement. This is derived form application of correspondence principle to elastic solution obtained by Sneddon 1965 (Int Journal of Engng Sci. vol3, pp.47-57), see eq 6.1
#MODEL PARAMETERS#


#CANTILEVER PARAMETERS#
A_free = params_simul[0]*1.0e-9  #free amplitude in nanometers
k1 = params_simul[1]
Q1 = params_simul[2]
freq = params_simul[10]*1.0e3
omega = 2*np.pi*freq
T = 1.0/freq
#CANTILEVER PARAMETERS#


#SIMULATION RESULT PARAMETERS
A = params_simul[8]*1.0e-9  #Reduced amplitude
Phase = params_simul[9]   #Output phase
nt = 10000 #number of steps for the total time of the analytical force solution
mT = 1  #number of periods used for the solution
#SIMULATION RESULT PARAMETERS


model_params = np.loadtxt('Model_params.txt', skiprows=1)
tau = model_params[:,0]
G = model_params[:,1]
Ge = model_params[0,2]




#RESULTS FROM THE SIMULATION#
simul = np.loadtxt('compu_40.00.txt', skiprows=1)
t_simul = simul[:,0]*1.0e-6
tip_simul = simul[:,1]*1.0e-9
F_simul = simul[:,2]*1.0e-9
t_simul[:] = t_simul[:]-t_simul[0]
#RESULTS FROM THE SIMULATION#


Zeq = z_eq(F_simul, tip_simul, A, 1)

t_an, t_in, t_fin_2, Tip_prescribed, F_bp, F_ss, F_tr, F_Gn, F_VdW = \
force_by_parts((Zeq), A, R, freq, G, tau, Ge, nu, H, nt)
#time shift of the sine with respect to the simulated tip#
t_shift = np.pi/(2*omega) - Phase*(np.pi/180.0)/omega
t_an[:] = t_an[:] - t_shift
#time shift of the sine with respect to the simulated tip#

F_Gn[F_Gn==0]=np.nan
F_tr[F_tr==0]=np.nan
F_ss[F_ss==0]=np.nan





#MAKING FIGURE
fig, ax = plt.subplots(2,1, figsize=(12,10))
figannotate(text='(a)', fs=20, pos=(-0.05,1.1))
plt.subplots_adjust(left=None, bottom=None, right=0.75, top=None, wspace=None, hspace=0.4)


# Twin the x-axis twice to make independent y-axes.
axes = [ax[0], ax[0].twinx(), ax[0].twinx()]
axes2  = [ax[1], ax[1].twinx()]

#getting an independent axis (NOT twin!)
axes3 = fig.add_axes(axes2[1].get_position(), frameon=False)

# Move the last y-axis spine over to the right by 15% of the width of the axes
axes[2].spines['right'].set_position(('axes', 1.15))

# To make the border of the right-most axis visible, we need to turn the frame
# on. This hides the other plots, however, so we need to turn its fill off.
axes[2].set_frame_on(True)
axes[2].patch.set_visible(False)



# Move the last y-axis spine over to the right by 20% of the width of the axes
axes3.spines['right'].set_position(('axes', 1.15))
#axes3.spines['top'].set_position(('axes', 0.9))

# To make the border of the right-most axis visible, we need to turn the frame
# on. This hides the other plots, however, so we need to turn its fill off.
axes3.set_frame_on(True)
axes3.patch.set_visible(False)



axes[0].plot(t_an*1.0e6, F_bp*1.0e9, color='0.6', lw=6, label='Force_analytical')
axes[0].plot(t_simul*1.0e6, F_simul*1.0e9,color='k', linestyle='solid', lw=2,  label='Force_simulation')
axes[0].set_ylabel("Total Force, nN", color='k', fontsize='16',fontweight='bold')
axes[0].set_xlabel('Time, a.u.', color='k', fontsize='16',fontweight='bold')
axes[0].legend(loc=3, fancybox=True, framealpha=0, fontsize=12)
axes[0].set_ylim(-20,60)
axes[0].set_xlim(0.4*T*1.0e6, 0.9*T*1.0e6)
axes[0].xaxis.set_ticks_position('bottom')

axes[1].plot(t_an*1.0e6, F_Gn*1.0e9, 'b--', lw=4, label='Relaxation Modulus')
axes[1].plot(t_an*1.0e6, (F_tr)*1.0e9, 'b', linestyle='dotted', lw=5, label='Transients')
axes[1].tick_params(axis='y', colors = 'blue')
axes[1].spines['right'].set_color('blue')
axes[1].spines['right'].set_linewidth(3.0)
axes[1].set_ylim(-10000,10000)
axes[1].set_ylabel("Transient and G(t), nN", color='b', fontsize='16',fontweight='bold')
axes[1].legend(loc=4, fancybox=True, framealpha=0, fontsize=12)
axes[1].xaxis.set_ticks_position('bottom')

axes[2].plot(t_an*1.0e6, (F_ss)*1.0e9, 'g', lw=4, ls='dashdot', label='Steady State')
axes[2].tick_params(axis='y', colors = 'green')
axes[2].spines['right'].set_color('green')
axes[2].spines['right'].set_linewidth(3.0)
axes[2].set_ylabel("Steady State, nN", color='g', fontsize='16',fontweight='bold')
axes[2].set_ylim(-1000,1000)
axes[2].legend(loc=2, fancybox=True, framealpha=0, fontsize=11)
axes[2].xaxis.set_ticks_position('bottom')


axes2[0].plot(Tip_prescribed*1.0e9, F_bp*1.0e9, color='0.6', lw=6, label='Force_analytical')
axes2[0].plot((tip_simul-Zr)*1.0e9, F_simul*1.0e9,color='k', linestyle='solid', lw=2,  label='Force_simulation')
axes2[0].set_ylabel("Total Force, nN", color='k', fontsize='16',fontweight='bold')
axes2[0].set_xlabel(r'$z_{t-s}, nm$', color='k', fontsize='20',fontweight='bold')
axes2[0].set_ylim(-20,60)
axes2[0].set_xlim(-4,10)
axes2[0].xaxis.set_ticks_position('bottom')

axes2[1].plot(Tip_prescribed*1.0e9, F_Gn*1.0e9, 'b--', lw=4, label='Relaxation Modulus')
axes2[1].plot(Tip_prescribed*1.0e9, (F_tr)*1.0e9, 'b', linestyle='dotted', lw=5, label='Transients')
axes2[1].tick_params(axis='y', colors = 'blue')
axes2[1].spines['right'].set_color('blue')
axes2[1].spines['right'].set_linewidth(3.0)
axes2[1].set_ylabel("Transient and G(t), nN", color='b', fontsize='16',fontweight='bold')
axes2[1].set_ylim(-10000,10000)
axes2[1].xaxis.set_ticks_position('bottom')


axes3.plot(Tip_prescribed*1.0e9, (F_ss)*1.0e9, 'g', lw=4, ls='dashdot', label='Steady State')
axes3.tick_params(axis='y', colors = 'green')
axes3.spines['right'].set_color('green')
axes3.spines['right'].set_linewidth(3.0)
axes3.tick_params(axis='x', colors = 'green')


axes3.spines['bottom'].set_linewidth(0.0)
axes3.spines['left'].set_linewidth(0.0)
axes3.set_yticks(np.arange(-1000, 1000, 500))
axes3.set_xlim(-20,300)
axes3.set_ylim(-1000,1000)
axes3.set_ylabel("Steady State, nN", color='g', fontsize='16',fontweight='bold')
axes3.yaxis.set_label_position('right')
axes3.yaxis.tick_right()
axes3.xaxis.tick_top()
axes3.spines['top'].set_color('green')
axes3.spines['top'].set_linewidth(3.0)



plt.savefig('F_bp.png', bbox_inches='tight')