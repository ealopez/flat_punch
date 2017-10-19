# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 13:31:54 2016

@author: Enrique Alejandro
"""

import numpy as np
from pylab import gcf

def G_lost(omega, G, tau, Ge = 0.0):
    "this function returns the value of G_loss for either a point value or an array of omega"
    if np.size(omega) == 1:  #calculation of point value of G_loss
        G_biprime = 0.0
        if np.size(G)>1: #the model has more than one arm
            G_biprime = sum(  G[:]*omega*tau[:]/( 1.0+pow(omega,2)*pow(tau[:],2) )  )
        else:  #the modela is the SLS
            G_biprime = (  G*omega*tau/( 1.0+pow(omega,2)*pow(tau,2) )  )
    else: #calculation of G_loss for an array of omega
        G_biprime = np.zeros(np.size(omega))
        if np.size(G) > 1: #the model has more than one arm
            for j in range(np.size(omega)):
                G_biprime[j] = sum(  G[:]*omega[j]*tau[:]/( 1.0+pow(omega[j],2)*pow(tau[:],2) )  )
        else: # the model is the SLS
            for j in range(np.size(omega)):
                G_biprime[j] = (  G*omega[j]*tau/( 1.0+pow(omega[j],2)*pow(tau,2) )  )              
    return G_biprime

def G_store(omega, G, tau, Ge = 0.0):
    "this function returns the value of G_store for either a point value or an array of omega"
    if np.size(omega) == 1:  #calculation of point value of G_loss
        G_prime = 0.0
        if np.size(G)>1: #the model has more than one arm
            Gg = Ge+sum(G[:])
            G_prime = Gg - sum(  G[:]/( 1.0+pow(omega,2)*pow(tau[:],2) )  )
        else:  #the modela is the SLS
            Gg = Ge+G
            G_prime = Gg - (  G/( 1.0+pow(omega,2)*pow(tau,2) )  )
    else: #calculation of G_loss for an array of omega
        G_prime = np.zeros(np.size(omega))
        if np.size(G) > 1: #the model has more than one arm
            Gg = Ge+sum(G[:])
            for j in range(np.size(omega)):
                G_prime[j] = Gg - sum(  G[:]/( 1.0+pow(omega[j],2)*pow(tau[:],2) )  )
        else: # the model is the SLS
            for j in range(np.size(omega)):
                Gg = Ge+G
                G_prime[j] = Gg - (  G/( 1.0+pow(omega[j],2)*pow(tau,2) )  )             
    return G_prime    


def G_abs(G_storage, G_loss):
    "returns the absolute modulus"
    G_absolute = np.sqrt(G_loss**2+G_storage**2)
    return G_absolute

def theta_loss_G(omega, G, tau, Ge = 0):
    "this function returns the loss angle from Generalized Maxwell Prony Coefficients"
    Gloss = G_lost(omega, G, tau, Ge)
    Gstorage = G_store(omega, G, tau, Ge)
    theta = np.arctan(Gloss/Gstorage)*180.0/np.pi
    return theta

def force_by_parts(Z, A, R, freq, G, tau, Ge = 0, nu=0.5, H=2.0e-19, nt=10000, mT=1.0):
    omega = 2.0*np.pi*freq
    alfa = 4.0*R/(1.0-nu)
    a = 0.2e-9   #intermolecular distance
    G_storage = G_store(omega, G, tau, Ge)
    G_loss = G_lost(omega, G, tau, Ge)
    G_absolute = G_abs(G_storage, G_loss)
    gama = H/(24.0*np.pi*a**2)   #surface energy
    F_VdW = - 4.0*np.pi*gama*R
    F_ss = np.zeros(nt)
    C = np.sqrt(A**2-Z**2)
    t_fin = 1.0/omega*(  np.arcsin(Ge*(Z/A)/G_absolute) - np.arctan(G_loss/G_storage) + 2*np.pi )
    t_fin_2 = t_fin
    theta = theta_loss_G(omega, G, tau, Ge)*np.pi/180.0
    N = np.size(G)  #number of Maxwell arms
    an = np.zeros(N)
    F_an = np.zeros(nt)
    F_Gn = np.zeros(nt)
    F_bp = np.zeros(nt)
    T = 2.0*np.pi/omega
    t = np.linspace(0,mT*T, nt)    #time array
    t_in = (np.arcsin(Z/A)+np.pi)/omega   #initial time at which the tip encounters the sample
    Tip_prescribed = np.zeros(nt)
    Tip_prescribed[:] = A*np.sin(omega*t[:]) + Z  #prescribed trajectory of the tip    
    for i in range (1,nt): #advancing in time
        sum_FGn = Ge
        sum_an = 0.0
        for j in range(0,N):
            an[j] = G[j]*(Z-C*omega*tau[j])/(1.0+omega**2*tau[j]**2)
            sum_an = sum_an + an[j]*np.exp(-(t[i]-t_in)/tau[j])
            sum_FGn = sum_FGn + G[j]*np.exp(-(t[i]-t_in)/tau[j])
        F_ss[i] = -alfa*A*G_absolute*np.sin(omega*t[i]+theta)
        F_an[i] = alfa*sum_an
        F_Gn[i] = -alfa*Z*(sum_FGn)
        if (t[i]>t_in):
            if (t[i]>t_fin_2):
                F_bp[i] = 0.0
            else:
                F_bp[i] = F_ss[i] + F_an[i]  + F_Gn[i] + F_VdW 
                if (F_bp[i] < F_VdW):  #SEMIANALYTICAL CALCULATION OF ACTUAL TIME WHERE THE ANALYTICAL FORCE CHANGES SIGN
                    t_fin_2 = t[i]  #Getting t double prime (the time when force becomes zero)
        else:
            F_bp[i] = 0.0
            F_an[i] = 0.0
            F_Gn[i] = 0.0            
    return t, t_in, t_fin_2, Tip_prescribed, F_bp, F_ss, F_an, F_Gn, F_VdW

def Ediss_by_parts(t_in, t_fin_2, Zeq, A, R, freq, G, tau, Ge = 0, nu=0.5, H=2.0e-19, nt=10000, mT=1.0):
    """this function decouples the contributions of dissipated energy from different sources: steady-state, transients, 
    relaxation modulus and van der Waals, as decoupled in analytical Eq for publication (Phys REvB)"""
    a = 0.2e-9 #interatomic distance
    omega = 2.0*np.pi*freq
    alfa = 4.0*R/(1.0-nu)
    Cn = np.sqrt(A**2-Zeq**2)
    Bn = A/Zeq*np.sin(omega*t_fin_2)
    N = np.size(G)
    an = np.zeros(N)
    an[:]= G[:]*(Zeq-Cn*omega*tau[:])/(1.0+tau[:]**2*omega**2)
    dn = np.zeros(N)
    Dn = np.sqrt(1.0-(Bn/A*Zeq)**2)
    dn[:]= tau[:]/(tau[:]**2*omega**2+1.0)*( np.exp((t_in-t_fin_2)/tau[:])*(tau[:]*omega*Bn*Zeq/A - Dn)+tau[:]*omega*Zeq/A - Cn/A)
    gamma = H/(24.0*a**2*np.pi)
    
    G_loss = G_lost(omega, G, tau, Ge)
    G_storage = G_store(omega, G, tau, Ge)
    E_loss_n = alfa*pow(A,2)/2.0*G_loss*Zeq/A*( A/Zeq*omega*(t_fin_2-t_in) - Cn/A + Bn*Dn )
    E_storage_n = alfa*pow(A,2)/2.0*G_storage*(Zeq/A)**2*(Bn**2-1)
    E_Ftransient_n = -A**2*alfa*omega/A*sum(an[:]*dn[:])
    E_F_relax_n = A**2*alfa*omega*Zeq/A*sum(G[:]*dn[:]) + alfa*A**2*(Zeq/A)**2*Ge*(Bn+1.0)  
    E_vdw_n = alfa*A**2* (4.0*np.pi*R*gamma/(alfa*Zeq))*(Zeq/A)**2*(Bn+1.0)
    E_total_n = E_loss_n + E_storage_n + E_Ftransient_n + E_F_relax_n + E_vdw_n
    return E_total_n, E_loss_n, E_storage_n, E_Ftransient_n, E_F_relax_n, E_vdw_n

def Vts_by_parts(t_in, t_fin_2, Zeq, A, R, freq, G, tau, Ge = 0, nu=0.5, H=2.0e-19, nt=10000, mT=1.0):
    """this function decouples the contributions of virial from different sources: steady-state, transients, 
    relaxation modulus and van der Waals, as decoupled in analytical Eq for publication (Phys REvB)"""
    N = np.size(G)
    a = 0.2e-9 #interatomic distance
    omega = 2.0*np.pi*freq
    alfa = 4.0*R/(1.0-nu)
    Cn = np.sqrt(A**2-Zeq**2)
    Bn = A/Zeq*np.sin(omega*t_fin_2)
    an = np.zeros(N)
    an[:]= G[:]*(Zeq-Cn*omega*tau[:])/(1.0+tau[:]**2*omega**2)
    dn = np.zeros(N)
    Dn = np.sqrt(1.0-(Bn/A*Zeq)**2)
    dn[:]= tau[:]/(tau[:]**2*omega**2+1.0)*( np.exp((t_in-t_fin_2)/tau[:])*(tau[:]*omega*Bn*Zeq/A - Dn)+tau[:]*omega*Zeq/A - Cn/A)
    gamma = H/(24.0*a**2*np.pi)
    gn = np.zeros(N)
    gn[:]= tau[:]/(tau[:]**2*omega**2+1.0)*( -np.exp((t_in-t_fin_2)/tau[:])*(tau[:]*omega*Dn + Bn*Zeq/A) - tau[:]*omega*Cn/A - Zeq/A )
    
    G_loss = G_lost(omega, G, tau, Ge)
    G_storage = G_store(omega, G, tau, Ge)
    vts_st = -alfa*A**2* G_storage/(4.0*np.pi)*Zeq/A*(A/Zeq*omega*(t_fin_2-t_in) + Cn/A - Bn*Dn)
    vts_loss = -alfa*A**2*G_loss/(4.0*np.pi)*(Zeq/A)**2*(Bn**2-1.0)
    vts_tr = alfa*A**2*omega/(2.0*np.pi*A)*sum(an[:]*gn[:])
    vts_Ge = alfa*A**2*Ge*Zeq/A*1.0/(2*np.pi)*(Cn/A+Dn) - alfa*A**2*omega/(2*np.pi)*Zeq/A*sum(G[:]*gn[:])
    vts_vdw = alfa*A**2*(2.0*R*gamma)/(A*alfa)*(Cn/A+Dn)
    vts = vts_st + vts_loss + vts_tr + vts_Ge + vts_vdw
    return vts, vts_loss, vts_st, vts_tr, vts_Ge, vts_vdw


            
def z_eq(force, z_ts, A, N):
    """this function receives three arrays: time, force and tip position. Those are results from tapping
    simulation over N periods"""
    nf = np.size(force)
    nnf = int(nf/N)  #getting only one period of the solution int(nf/number of periods in solution)
    F_sim = force[:nnf]
    tip_sim = z_ts[:nnf]
    flag = 0
    b=0
    c=0
    for i in range(1,nnf-1):
        if (F_sim[i] < F_sim[i-1] and F_sim[i] < F_sim[i+1] and flag==0):
            b = i
            flag = 1
    flag =0
    for i in range(1,nnf-1):
        if (tip_sim[i] < tip_sim[i-1] and tip_sim[i] < tip_sim[i+1] and flag==0):
            c = i
            flag = 1
    indent = tip_sim[b] - tip_sim[c]
    if (b==0):
        Zeq = np.nan
    else:
        Zeq = A-indent  
    return Zeq         
        
def periodic_extension(t, f_t, N=25):
    """this function receives one single period of a periodic function and returns
    a periodic extension over a defined number of periods, the default is 25"""
    period = t[np.size(t)-1]
    #dt = time[1]-time[0]
    t_extend = np.linspace(0,N*period, np.size(t)*N)  #extension of time array
    f_t_extend = np.zeros(np.size(t_extend))
    for i in range(np.size(f_t_extend)):
        if i < np.size(t):
            f_t_extend[i] = f_t[i]
        else:
            f_t_extend[i] = f_t_extend[i-np.size(t)]         
    return t_extend, f_t_extend
    
def G_t(t, G, tau, Ge = 0.0):
    """this function returns the relaxation modulus in time"""
    G_rel = np.zeros(np.size(t))
    if np.size(G) == 1:  #the model is the SLS
        for i in range(np.size(t)):
            G_rel[i] = Ge + G*np.exp(-t[i]/tau)
    else: #the model has more than one arm
        for i in range(np.size(t)):
            G_rel[i] = Ge + sum(G[:]*np.exp(-t[i]/tau[:]))
    return G_rel

    

def Amp_Phase(t, f_t, freq):
    """this function calculates amplitude and phase using the in-phase and in-quadrature integrals"""
    if t[0] > 0.0:
        t[:] = t[:t]-t[0]
    dt = t[1] - t[0]
    I = 0.0
    K = 0.0
    for i in range(np.size(f_t)):
        I = I + f_t[i]*np.cos(2.0*np.pi*freq*t[i])*dt
        K = K + f_t[i]*np.sin(2.0*np.pi*freq*t[i])*dt
    Amp = 1.0/(t[np.size(t)-1])*np.sqrt(I**2+K**2) *2.0
    Phase = np.arctan(K/I)*180.0/np.pi
    return Amp, Phase


    
def E_diss(z, Fts):
    """Output: this function calculates the tip-sample dissipation"""
    """Input: tip position and tip-sample force arrays"""
    Ediss = 0.0
    for i in range(1,len(z)-1):
        Ediss = Ediss - Fts[i]*(z[i+1]-z[i-1])/2.0   #based on integral of Fts*dz/dt*dt, dz/dt=(z[i+1]-z[i-1])/(2.0*dt) Central difference approx
    return Ediss

def V_ts(z, Fts, dt):
    """Output: virial"""
    """Input: tip position and tip-sample force arrays, and timestep"""
    Vts = 0.0
    for i in range(len(z)):
        Vts = Vts + Fts[i]*z[i]*dt
    return Vts/(dt*len(z))     #virial is 1/T*S(Fts*z*dt) from 0 to T, being T total experimental time

def av_Fts(Fts, dt):
    avF = 0.0
    for i in range(len(Fts)):
        avF = avF + Fts[i]*dt
    return avF/(dt*len(Fts))



def figannotate(text='a',fs=8,ff='helvetica',fw='bold',pos=(-0.02,1),ax=0,ha='right',va='top'):
    if ax==0:
        ax = gcf().axes
    for i in range(len(ax)):
        ax[i].text(pos[0],pos[1],text,ha=ha,va=va,fontsize=fs,family=ff,weight=fw,transform=ax[i].transAxes)
        if len(text)>2:
            text = text[0]+'%s'%(chr(ord(text[1])+1))+text[2:]
        else:
            text = '%s'%(chr(ord(text[0])+1))+text[1:]
