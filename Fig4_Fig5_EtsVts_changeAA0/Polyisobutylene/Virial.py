# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:36:01 2017

@author: enrique
test: getting Zeq to get analytical Ediss and comparing to numerical Ediss
"""

from flatpunch import z_eq, force_by_parts, Vts_by_parts
import numpy as np
import matplotlib.pyplot as plt
from glob import glob


#OPENING GENERAL FILES
summ = np.loadtxt('summary.txt', skiprows = 1)
Amp_free = summ[0,0]*1.0e-9
Amp = summ[:,8]*1.0e-9  #reduced tapping amplitude in meters
A_ratio = Amp/Amp_free
zq = summ[:,6]   #nominal zeq from the beggining of simulation
zlist = zq.tolist()  #transfering nominal Zeq to list, so it can be compared using index function
zq = zq.tolist()
zr = summ[:,7]
zs = (zq-zr)  #z equilibrium according to the c code
indent = (Amp - zs)  #total indentation
R = summ[0,3]*1.0e-9 #tip radius in m
freq = summ[0,10]*1000.0  #tapping frequency in Hz
nu = summ[0,4]
Virial_compu = summ[:,12]*1.0e-18  #dissipated energy calculated in the c code (numerically)

material = np.loadtxt('Model_params.txt', skiprows=1)
tau = material[:,0]
G = material[:,1]
Ge = material[0,2]
#OPENING GENERAL FILES

t_in =0.0
t_fin = 0.0

files = glob('compu_*.txt')  #open all files that have that extension
Vts_an = np.zeros(np.size(files))
Vts_loss= np.zeros(np.size(files))
Vts_store= np.zeros(np.size(files))
Vts_tr= np.zeros(np.size(files))
Vts_relax= np.zeros(np.size(files))
Vts_vdw= np.zeros(np.size(files))
zeq = np.zeros(np.size(files))

i=0
for f in files:
    sim = np.loadtxt(files[i], skiprows=1)
    tip = sim[:,1]*1.0e-9
    force = sim[:,2]*1.0e-9
    for j in files[i]:  #run through file's name to get numebr in it
        b = files[i].index('x')
        c = files[i].index('_')
    number= int(files[i][c+1:b-5]) #getting number in the file's name
    k = zlist.index(number) #position where nominal equilibrium position coincides with
    zeq[k] = z_eq(force, tip, Amp[k], 25.0)
    _, t_in, t_fin, _, _, _, _, _, _ = \
    force_by_parts(zeq[k], Amp[k], R, freq, G, tau, Ge, nu, 2.0e-19, 10000, 1.0)
    Vts_an[k], Vts_loss[k], Vts_store[k], Vts_tr[k], Vts_relax[k], Vts_vdw[k] = Vts_by_parts(t_in, t_fin, zeq[k], Amp[k], R, freq, G, tau, Ge, nu, 2.0e-19, 10000, 1.0)
    i=i+1



#plotting figure:
fig, ax = plt.subplots(figsize=(10.1,5))


# Twin the x-axis twice to make independent y-axes.
axes = [ax, ax.twinx(), ax.twinx()]

# Make some space on the right side for the extra y-axis.
fig.subplots_adjust(right=0.75)

# Move the last y-axis spine over to the right by 20% of the width of the axes
axes[-1].spines['right'].set_position(('axes', 1.15))

# To make the border of the right-most axis visible, we need to turn the frame
# on. This hides the other plots, however, so we need to turn its fill off.
axes[-1].set_frame_on(True)
axes[-1].patch.set_visible(False)

# And finally we get to plot things...
axes[0].plot(A_ratio, Virial_compu*1.0e18,'-^',color='dimgray', markersize=8.0, lw=1.0, label = 'Simulation')
axes[0].plot(A_ratio, Vts_an*1.0e18, 'k-*', markersize=10.0, lw=1.0, label = 'Analytical')
axes[0].set_ylabel('Total virial, aJ', color='k', fontsize='18')
axes[0].tick_params(axis='y', colors = 'k')
axes[0].spines['left'].set_color('k')
axes[0].spines['left'].set_linewidth(3.0)
axes[0].set_ylim(-200,20)
axes[0].legend(loc=2, fancybox=True, framealpha=0, fontsize=11)


axes[1].set_ylabel('Steady State, aJ', color='blue', fontsize='18')
axes[1].plot(A_ratio, Vts_loss*1.0e18, 'b-o',markersize=7.0, lw=1.0, label = 'Loss')
axes[1].plot(A_ratio, Vts_store*1.0e18, 'b-D',markersize=7.0, lw=1.0, label = 'Storage')
axes[1].tick_params(axis='y', colors = 'blue')
axes[1].spines['right'].set_color('blue')
axes[1].spines['right'].set_linewidth(3.0)
axes[1].set_ylim(-3000,500)
axes[1].legend(loc=3, fancybox=True, framealpha=0, fontsize=11)


axes[2].plot(A_ratio, Vts_tr*1.0e18, 'g-x', markeredgewidth=2.0, markersize=8.0, lw=1.0, label = 'Transient')
axes[2].plot(A_ratio, Vts_relax*1.0e18, 'g-s', lw=1.0,markersize=7.0, label = 'Equilibrium')
axes[2].set_ylabel('Transient and Relaxation, aJ', color='green', fontsize='18')
axes[2].tick_params(axis='y', colors = 'green')
axes[2].spines['right'].set_color('green')
axes[2].spines['right'].set_linewidth(3.0)
axes[2].set_ylim(-10000,10000)
axes[2].legend(loc=4, fancybox=True, framealpha=0, fontsize=11)


axes[0].set_xlabel(r'$A/A_0$', color='k', fontsize='20',fontweight='bold')
axes[0].set_xlim(0.1,1.0)


plt.tight_layout(pad=0.4, w_pad=2.0, h_pad=0.5)
plt.savefig('Virial.png', bbox_inches='tight')