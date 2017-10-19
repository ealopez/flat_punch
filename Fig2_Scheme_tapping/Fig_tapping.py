# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 09:37:42 2016

@author: Enrique Alejandro
"""

import numpy as np
import matplotlib.pyplot as plt

amph = np.loadtxt('summary.txt', skiprows=1)
amp = amph[3]*10.0*1.0e-9 #multiplying by free amplitude and then converting to m
phi = amph[4]*np.pi/180.0  #converting phase to radians
omega = 2*np.pi*6.0e5
compu = np.loadtxt('compu_4.00.txt', skiprows=1)
t = compu[:,0]*1.0e-6  #converting to seconds
tip = compu[:,1]*1.0e-9  #converting to meters
cos_ref = compu[:,2]*1.0e-9  #converting to meters
Fts = compu[:,3]*1.0e-9  #converting to Newtons
xb = compu[:,4]*1.0e-9 #converting to meters
Zeq = np.zeros(np.size(tip))
Zeq[:] = 4.0e-9



fig, ax = plt.subplots(1,1,figsize=(12,5))
ax.plot(t*1.0e6, tip*1.0e9, 'b', ls='dashdot', lw=2, label=r'$z_{t-s}(t)=Z_{eq}+z(t)$')
ax.plot(t*1.0e6, xb*1.0e9, 'g', lw=3, label=r'$Sample Position$')
ax.legend(loc=4, fontsize=18, frameon=False)
ax.plot(t*1.0e6, cos_ref*1.0e9, 'k', lw=1, ls='dotted', label=R'$z(t)=A*sin(\omega*t)$')
ax.plot(t*1.0e6, Zeq*1.0e9, color='dimgray', ls='dashed', lw=1, label=r'$Z_{eq}$')
ax.legend(loc=4, fontsize=16, frameon=False)
#plt.plot(t*1.0e6, amp*np.cos(omega*t-phi)*1.0e9, 'c', lw=2, label='Sine Reference')

plt.xlabel('time, a.u.', fontsize='20',fontweight='bold')
plt.ylabel('Z-position, nm',fontsize='20',fontweight='bold')
ax.set_xlim(840.2,839.95+(2*np.pi/omega)*1.0e6*1.4)
ax.set_ylim(-12,13.5)
ax.set_xticks([])
#ax.set_yticks([])
#ax.axis('off')

plt.savefig('Tapping_Scheme.png', bbox_inches='tight')