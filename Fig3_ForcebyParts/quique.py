# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 13:31:54 2016

@author: Enrique Alejandro
"""

from pylab import*
from aNoh import*
import numpy as np
from lmfit import minimize, Parameters

def smear(x, t, tr=0.1, st=1.0):
    "this function smears an array providing as reference also time array, time resolution and total time"
    nt = size(t)
    prints = 0
    i =0
    x_smear = []
    t_smear = []
    while i < (nt):
        if t[i] >= prints*tr and t[i]<=(st+tr) :
            x_smear.append(x[i])
            t_smear.append(t[i])
            prints = prints + 1
        i = i + 1
    return asarray(x_smear), asarray(t_smear)

a = linspace(3.0,6.0,10.0e3)
tiempo = linspace(0,1.0,10000)
a_smear, tiempo_smear = smear(a, tiempo, 1.0e-3, 1.0)

def log_tw(de0, maxi, nn=10):
    "this function generates a frequency or time array in log scale"
    epsilon = []
    w = de0
    de = de0
    #epsilon[0] = 0
    prints = 1
    epsilon.append(de0)
    while w < maxi:
        w = w + de
        if w < maxi:
            epsilon.append(w)              
        prints = prints + 1 
        if prints == nn:
            de = de*10
            prints = 1    
    return asarray(epsilon)


def log_scale(x, t, tr=0.1, st=1.0, nn = 10):
    "this receives an array and write it weighting it in logarithmic scale"
    prints = 1
    nt = size(x)
    i =0
    x_log = []
    t_log = []
    while i <nt:
        if t[i] >= prints*tr and t[i]<=st :
            x_log.append(x[i])
            t_log.append(t[i])
            prints = prints + 1
        i = i + 1
        if prints == nn:
            tr = tr*10
            prints = 1
    return asarray(x_log), asarray(t_log)
    
def linear_force(params, t, force):
    "this function returns the"
    p = params.valuesdict()
    A = p['A']
    model = (A)*t
    return  (model - force) #calculating the residual

def logdata(x):
    "this functions takes the array x and weight them in logarithmic scale"
    nt = size(x)
    x_r = []
    counter = 1
    print_counter =0
    i = 0
    gap = 1
    #flag = 0
    #prints = 0
    while i<nt:
        i = i + gap    
        gap = (10**counter)/10.0   #space between prints
        print_counter = print_counter+1
        if print_counter == 9:
            print_counter = 0
            counter = counter + 1
        if (i < (nt)) :
            #if (flag==0):  #this make sure the first point of the interesting range is printed (smaller timestep given by time resolution)
                #flag = 1
                #x_r.append(x[i-gap])
                #prints = prints +1
            x_r.append(x[i])
            #prints = prints + 1
    y = asarray(x_r)
    return  y

        
def percent_error(th, ex):
    "this functions calculates the percentual error"     
    error = abs((th-ex)/th)*100.0
    return error

def abs_error(th, ex):
    "this functions calculates the percentual error"     
    error = abs((th-ex))
    return error


#CALCULATION OF HARMONIC QUANTITIES#
def J_storage(omega, Jg, J, tau):
    "this function gives an array of J_storage on omega"
    J_prime = zeros(size(omega))
    for i in range(size(omega)):
        if size(J) > 1:
            J_prime[i] = Jg + sum( J[:]/(1.0 + (  pow(omega[i],2)*pow(tau[:],2) ) ) )
        else:
            J_prime[i] = Jg + ( J/(1.0 + (  pow(omega[i],2)*pow(tau,2) ) ) )
            
    return J_prime

def J_loss(omega, Jg, J, tau, phi = 0):
    "this function gives an array of J_loss on omega"
    J_biprime = zeros(size(omega))
    for i in range(size(omega)):
        if size(J)>1:
            J_biprime[i] = sum( J[:]*omega[i]*tau[:]/(1.0 + (pow(omega[i],2)*pow(tau[:],2)) ) ) + phi/omega[i]
        else:
            J_biprime[i] = ( J*omega[i]*tau/(1.0 + (pow(omega[i],2)*pow(tau,2)) ) ) + phi/omega[i]
    return J_biprime

def theta_loss(omega, Jg, J, tau, phi =0):    
    "this function returns the loss angle from Generalized Kelvin Voigt Prony Coefficients"
    Jloss = J_loss(omega, Jg, J, tau, phi =0)
    Jstorage =  J_storage(omega, Jg, J, tau)
    theta = arctan(Jloss/Jstorage)*180/pi
    return theta

def theta_loss_G(omega, G, tau, Ge = 0):
    "this function returns the loss angle from Generalized Maxwell Prony Coefficients"
    Gloss = G_lost(omega, G, tau, Ge)
    Gstorage = G_store(omega, G, tau, Ge)
    theta = arctan(Gloss/Gstorage)*180.0/pi
    return theta

def chi_th(t, Jg, J, tau, phi = 0):
    "this function gives the strain response to a unit slope stress"
    if (size(J)) > 1:
        Je = sum(J[:])+Jg
    else:
        Je = J+Jg
    chi = zeros(size(t))
    for i in range (size(t)):
        if (size(J)) >1 :
            chi[i] = Je*t[i] + sum(J[:]*tau[:]*(exp(-t[i]/tau[:])-1.0)) + 1/2.0*phi*pow(t[i],2)
        else:
            chi[i] = Je*t[i] + (J*tau*(exp(-t[i]/tau)-1.0)) + 1/2.0*phi*pow(t[i],2)
    return chi

def chi_exp(b, Fdot, h):
    return b/Fdot*pow(h, 1.5)
    
def compliance(t, Jg, J, tau, phi=0):
    "this function returns the compliance in time t, for a model with given set of parameters"
    if (size(J)) > 1:
        Je = sum(J[:])+Jg
    else:
        Je = J+Jg
    comp = zeros(size(t))
    for i in range (size(t)):
        if size(J)>1:
            comp[i] = Je - sum(J[:]*exp(-t[i]/tau[:])) + phi*t[i]
        else:
            comp[i] = Je - (J*exp(-t[i]/tau)) + phi*t[i]
    return comp

def J_loss_farti(Jstore_arti, Jloss_arti, E, b):
    "this function is designed for the flat-end indenter, it calculates the loss compliance from the artificial ones"
    return Jloss_arti/ ( pow((1.0-E/b*Jstore_arti),2) + pow((E/b*Jloss_arti),2) )

def J_storage_farti(Jstore_arti, Jloss_arti, E, b):
    "this function is designed for the flat-end indenter, it calculates the loss compliance from the artificial ones"
    return ( Jstore_arti - E/b*pow(Jstore_arti,2) - E/b*pow(Jloss_arti,2) ) / ( pow((1.0-E/b*Jstore_arti),2) + pow((E/b*Jloss_arti),2) )
   
def Ut(time, J, tau, phi=0):
    "this function gives the response of a unit strain impulse"
    U = zeros(size(time))
    for i in range(size(time)):
        U[i] = sum(J[:]/tau[:]*exp(-time[i]/tau[:])) + phi
    return U

def conv_U_F(time, F, Jg, J, tau, phi=0):
    "this function convolves force and U(t)"
    dt = av_dt(time)
    U = zeros(size(time))
    for i in range(size(time)):
        U[i] = sum(J[:]/tau[:]*exp(-time[i]/tau[:])) +phi
    conv = convolve(U, F, mode='full')*dt
    conv = conv[range(size(F))] + Jg*F
    return conv

def G_lost(omega, G, tau, Ge = 0.0):
    "this function returns the value of G_loss for either a point value or an array of omega"
    if np.size(omega) == 1:  #calculation of point value of G_loss
        G_biprime = 0.0
        if size(G)>1: #the model has more than one arm
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
        if size(G)>1: #the model has more than one arm
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
    G_absolute = sqrt(G_loss**2+G_storage**2)
    return G_absolute

def force_by_parts(Z, A, R, freq, G, tau, Ge = 0, nu=0.5, H=2.0e-19, nt=10000, mT=1.0):
    omega = 2.0*pi*freq
    alfa = 4.0*R/(1.0-nu)
    a = 0.2e-9   #intermolecular distance
    G_storage = G_store(omega, G, tau, Ge)
    G_loss = G_lost(omega, G, tau, Ge)
    G_absolute = G_abs(G_storage, G_loss)
    gama = H/(24.0*np.pi*a**2)   #surface energy
    F_VdW = - 4.0*np.pi*gama*R
    F_ss = zeros(nt)
    C = sqrt(A**2-Z**2)
    t_fin = 1.0/omega*(  arcsin(Ge*(Z/A)/G_absolute) - arctan(G_loss/G_storage) + 2*pi )
    t_fin_2 = t_fin
    theta = theta_loss_G(omega, G, tau, Ge)*np.pi/180.0
    N = size(G)  #number of Maxwell arms
    an = zeros(N)
    F_an = zeros(nt)
    F_Gn = zeros(nt)
    F_bp = zeros(nt)
    T = 2.0*pi/omega
    t = linspace(0,mT*T, nt)    #time array
    t_in = (arcsin(Z/A)+pi)/omega   #initial time at which the tip encounters the sample
    Tip_prescribed = zeros(nt)
    Tip_prescribed[:] = A*sin(omega*t[:]) + Z  #prescribed trajectory of the tip    
    for i in range (1,nt): #advancing in time
        sum_FGn = Ge
        sum_an = 0.0
        for j in range(0,N):
            an[j] = G[j]*(Z-C*omega*tau[j])/(1.0+omega**2*tau[j]**2)
            sum_an = sum_an + an[j]*exp(-(t[i]-t_in)/tau[j])
            sum_FGn = sum_FGn + G[j]*exp(-(t[i]-t_in)/tau[j])
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
    omega = 2.0*pi*freq
    alfa = 4.0*R/(1.0-nu)
    Cn = sqrt(A**2-Zeq**2)
    Bn = A/Zeq*sin(omega*t_fin_2)
    N = size(G)
    an = zeros(N)
    an[:]= G[:]*(Zeq-Cn*omega*tau[:])/(1.0+tau[:]**2*omega**2)
    dn = zeros(N)
    Dn = sqrt(1.0-(Bn/A*Zeq)**2)
    dn[:]= tau[:]/(tau[:]**2*omega**2+1.0)*( exp((t_in-t_fin_2)/tau[:])*(tau[:]*omega*Bn*Zeq/A - Dn)+tau[:]*omega*Zeq/A - Cn/A)
    gamma = H/(24.0*a**2*pi)
    
    G_loss = G_lost(omega, G, tau, Ge)
    G_storage = G_store(omega, G, tau, Ge)
    E_loss_n = alfa*pow(A,2)/2.0*G_loss*Zeq/A*( A/Zeq*omega*(t_fin_2-t_in) - Cn/A + Bn*Dn )
    E_storage_n = alfa*pow(A,2)/2.0*G_storage*(Zeq/A)**2*(Bn**2-1)
    E_Ftransient_n = -A**2*alfa*omega/A*sum(an[:]*dn[:])
    E_F_relax_n = A**2*alfa*omega*Zeq/A*sum(G[:]*dn[:]) + alfa*A**2*(Zeq/A)**2*Ge*(Bn+1.0)  
    E_vdw_n = alfa*A**2* (4.0*pi*R*gamma/(alfa*Zeq))*(Zeq/A)**2*(Bn+1.0)
    E_total_n = E_loss_n + E_storage_n + E_Ftransient_n + E_F_relax_n + E_vdw_n
    return E_total_n, E_loss_n, E_storage_n, E_Ftransient_n, E_F_relax_n, E_vdw_n

def Vts_by_parts(t_in, t_fin_2, Zeq, A, R, freq, G, tau, Ge = 0, nu=0.5, H=2.0e-19, nt=10000, mT=1.0):
    """this function decouples the contributions of virial from different sources: steady-state, transients, 
    relaxation modulus and van der Waals, as decoupled in analytical Eq for publication (Phys REvB)"""
    N = size(G)
    a = 0.2e-9 #interatomic distance
    omega = 2.0*pi*freq
    alfa = 4.0*R/(1.0-nu)
    Cn = sqrt(A**2-Zeq**2)
    Bn = A/Zeq*sin(omega*t_fin_2)
    an = zeros(N)
    an[:]= G[:]*(Zeq-Cn*omega*tau[:])/(1.0+tau[:]**2*omega**2)
    dn = zeros(N)
    Dn = sqrt(1.0-(Bn/A*Zeq)**2)
    dn[:]= tau[:]/(tau[:]**2*omega**2+1.0)*( exp((t_in-t_fin_2)/tau[:])*(tau[:]*omega*Bn*Zeq/A - Dn)+tau[:]*omega*Zeq/A - Cn/A)
    gamma = H/(24.0*a**2*pi)
    gn = zeros(N)
    gn[:]= tau[:]/(tau[:]**2*omega**2+1.0)*( -exp((t_in-t_fin_2)/tau[:])*(tau[:]*omega*Dn + Bn*Zeq/A) - tau[:]*omega*Cn/A - Zeq/A )
    
    G_loss = G_lost(omega, G, tau, Ge)
    G_storage = G_store(omega, G, tau, Ge)
    vts_st = -alfa*A**2* G_storage/(4.0*pi)*Zeq/A*(A/Zeq*omega*(t_fin_2-t_in) + Cn/A - Bn*Dn)
    vts_loss = -alfa*A**2*G_loss/(4.0*pi)*(Zeq/A)**2*(Bn**2-1.0)
    vts_tr = alfa*A**2*omega/(2.0*pi*A)*sum(an[:]*gn[:])
    vts_Ge = alfa*A**2*Ge*Zeq/A*1.0/(2*pi)*(Cn/A+Dn) - alfa*A**2*omega/(2*pi)*Zeq/A*sum(G[:]*gn[:])
    vts_vdw = alfa*A**2*(2.0*R*gamma)/(A*alfa)*(Cn/A+Dn)
    vts = vts_st + vts_loss + vts_tr + vts_Ge + vts_vdw
    return vts, vts_loss, vts_st, vts_tr, vts_Ge, vts_vdw



def Ediss_Tamayo(k, Q, A_free, A, Phase):
    Ediss = (pi*k*A**2/Q)*( (A_free/A)*sin(Phase*pi/180.0) - 1.0 )
    return Ediss

def virial_Lozano(k, Q, A_free, A, Phase):
    Vts = -(k*A*A_free)/(2.0*Q)*cos(Phase*pi/180.0)
    return Vts
    
def force_simp(mT, indentation, A, R, alfa, freq, G, tau, Ge=0, nt=10000.0, H=2.0e-19):
    "this function returns the analytical force of a single tapping flat-end punch interacting with a generalized Maxwell viscoelastic sample"
    "it contains the force in time as defined in: Derivation of Analytical Relations of an AFM Probe Tapping a Viscoelastic Surface Containing Multiple Characteristic Times"    
    "the analytical equation corresponds to the algebraically simplified version"    
    G_storage = G_store(freq*2*pi, G, tau, Ge)
    G_loss = G_lost(freq*2*pi, G, tau, Ge)
    a = 0.2e-9
    Zs = -(A-indentation)  #negative distance between the equilibrium position and the sample
    F_VdW = - H*R/(6.0*a**2)
    F_Ge = alfa*Zs*Ge
    T = 1.0/freq
    t = linspace(0,mT*T, nt)
    omega = 2*pi*freq
    t_in = (arcsin(-Zs/A)+pi)/omega
    Tip_prescribed = zeros(nt)
    Tip_prescribed[:] = A*sin(omega*t[:]) - Zs  #prescribed trajectory of the tip
    F_simp = zeros(nt)
    F_storage = zeros(nt)
    F_loss = zeros(nt)
    F_tr = zeros(nt)
    C = sqrt(1.0-pow(Zs/A,2))
    N = size(G)
    G_absolute = sqrt(G_loss**2+G_storage**2)
    t_fin = 1.0/omega*(  arcsin(Ge*(Zs/A)/G_absolute) - arctan(G_loss/G_storage) + 2*pi )
    t_fin_2 = t_fin
    sum_tr=0.0
    for i in range (1,nt): #advancing in time
        sum_tr = 0.0
        for j in range(0,N):
            sum_tr = sum_tr + (G[j]*tau[j])/(1.0 + (omega**2*tau[j]**2))*exp(-(t[i]-t_in)/tau[j])*(omega*tau[j]*Zs-A*C)
        F_storage[i] = alfa*(-G_storage*A*sin(omega*t[i]))
        F_loss[i] = alfa*(-G_loss*A*cos(omega*t[i]))
        F_tr[i] = alfa*omega*sum_tr
        if (t[i]>t_in):
            if (t[i]>t_fin_2):
                F_simp[i] = 0.0
            else:
                F_simp[i] = F_storage[i]+ F_loss[i] + F_tr[i] + F_VdW + F_Ge
                if (F_simp[i] < F_VdW):  #SEMIANALYTICAL CALCULATION OF ACTUAL TIME WHERE THE ANALYTICAL FORCE CHANGES SIGN
                    t_fin_2 = t[i]  #Getting t double prime (the time when force becomes zero)
        else:
            F_simp[i] = 0.0
    return t, F_simp

def sign_change(array):
    """this function gives the position at which the numbers in the array change sign for first time"""
    i = 0
    if array[0]<0:
        while (array[i]) < 0:
            i = i+1
    else:
        while (array[i]) > 0:
            i = i+1
    return i

def av_dt(array):
    "this function returns the average of the timesteps in a time array"
    i = 0
    k = 0.0
    for i in range(np.size(array)-1):
        k = k + (array[i+1]-array[i])
        dt = k/(np.size(array)-1)
    return dt
            
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
        Zeq = nan
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

def A_phi(t, f_t, h=1.0, f1=1):
    """This function gets the amplitude and phase Fourier Components of a defined harmonic
    from a function given in time"""
    if f1 ==1:  #if frequency of 1st harmonic is not provided
        f1=1.0/t[np.size(t)-1]
    if t[0] > 0.0:
        t[:] = t[:t]-t[0]
    f_nu = np.fft.fft(f_t)/(np.size(t))  #this transforms the time dependent function to its transformed frequency dependent
    nu = np.fft.fftfreq(np.size(t), t[1]-t[0])  #this takes care of finding the frequencies associated to the transformed function
    A = np.real(f_nu[0:np.size(t)/2.0+1])  #this the real part of the complex Fourier coefficients, only for the positive frequencies
    B = -np.imag(f_nu[0:np.size(t)/2.0+1])  #this the imaginary part of the complex Fourier coefficients, only for the positive frequencies
    nu_real = nu[0:np.size(t)/2.0+1]
    Amp = np.sqrt(A**2+B**2) *2  #Getting the amplitude of the Fourier Coefficient and multiplying it for two (the amplitude of the sinudoid is two time the amplitude of the Fourier Coefficient)
    Phase = np.arctan (B/A)*180.0/np.pi
    dnu = nu_real[1]-nu_real[0]
    index = int(f1/dnu)*h+h
    return Amp[index], Phase[index]
    

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


#NEED TO DEBUG:
def A_phi_array(t, f_t, f1=1):
    """This function gets all the available amplitude and phase Fourier Components of function given in time"""
    if f1 ==1:  #if frequency of 1st harmonic is not provided
        f1=1.0/t[np.size(t)-1]
    if t[0] > 0.0:
        t[:] = t[:t]-t[0]
    f_nu = np.fft.fft(f_t)/(np.size(t))  #this transforms the time dependent function to its transformed frequency dependent
    nu = np.fft.fftfreq(np.size(t), t[1]-t[0])  #this takes care of finding the frequencies associated to the transformed function
    A = np.real(f_nu[0:np.size(t)/2.0+1])  #this the real part of the complex Fourier coefficients, only for the positive frequencies
    B = -np.imag(f_nu[0:np.size(t)/2.0+1])  #this the imaginary part of the complex Fourier coefficients, only for the positive frequencies
    nu_real = nu[0:np.size(t)/2.0+1]
    dnu = nu_real[1]-nu_real[0]
    n = int(f1/dnu)
    A = A[::n+1]*1.0
    B = B[::n+1]*1.0
    freq = nu_real[::n+1]*1.0
    return A, B, freq

def f_harmonic(A, B, freq_0, t):
    """This function builds a function in time using trigonometric Fourier series with given coefficients"""
    for i in range(np.size(t)):
        suma = 0.0
        for j in range(1, np.size(A)):
            suma = suma + A[j]*np.cos(2.0*np.pi*j*freq_0*t[i]) +  B[j]*np.sin(2.0*np.pi*j*freq_0*t[i])
        force_t[i] = A[0]/2.0 + suma
    return force_t

def derivative_cd(f_t, t):
    """this function calculates the derivative of a given array using central difference scheme"""
    f_prime = np.zeros(np.size(f_t))
    for i in range(np.size(f_t)):  #calculation of derivative using central difference scheme
        if i == 0:
            f_prime[i] = (f_t[1]-f_t[0])/(t[1]-t[0])
        else:
            if i == np.size(f_t)-1:
                f_prime[np.size(f_t)-1] = (f_t[np.size(f_t)-1]-f_t[np.size(f_t)-2])/(t[np.size(f_t)-1]-t[np.size(f_t)-2])
            else:
                f_prime[i] = (f_t[i+1]-f_t[i-1])/(t[i+1]-t[i-1])
    return f_prime

def repulsive_FD(raw):
    """this function generates time, tip, and deflection array for the repulsive part only (force in time with positive slope)"""
    """Input: raw average data from Hanaul's averaging FD curves code, specifically raw = np.loadtxt('raw.txt') """
    """Output: file with time, tip and deflection arrays with only the repulsive part, also gives the deflection offset"""
    """Deflection offset is the lowest deflection right before force becomes repulsive (minimimum in the force in time curve) """
    t = raw[:,0]   #pulling time array
    tip = raw[:,2]  #pulling tip position
    D = raw[:,3]   #pulling deflection
    index_min = np.argmin(D)   #getting turning point of force in time
    defl = D[range(index_min, np.size(D))]
    def_offset = defl[0]
    defl = defl - def_offset
    tiempo = t[range(index_min, np.size(D))]
    avdt = av_dt(tiempo)
    tiempo = tiempo-tiempo[0] + avdt
    tip = tip[range(index_min, np.size(D))]
    tip = tip - tip[0]
    return tiempo, tip, defl, def_offset

def repulsive_array(indentation, defl, k=1):
    """This function returns the repulsive portion of an FD curve"""
    """Input: deflection, tip position arrays"""
    """Output: Tip Position and Force """
    index_min = np.argmin(defl) #getting turning point of force in time
    D = defl[range(index_min, np.size(defl))]
    def_offset = D[0]
    D = D - def_offset
    F = D*k
    tip = indentation[range(index_min, np.size(defl))]
    tip = tip - tip[0]
    return tip, F
    
def FD_log(x, f_x, nn=20, liminf = 1, suplim = 1):
    """this function returns time and f_time arrays equally spaced in logarithmic scale"""
    """Input: time starting with average dt(comming from repulsive_FD function), and f_time related to that time array"""
    if liminf ==1:
        lim_inf = round(np.log10(x[0]),2)
    else:
        lim_inf = round(np.log10(liminf),2)
    if suplim ==1:
        sup_lim = round(np.log10(x[np.size(x)-1]),2)
    else:
        sup_lim = round(np.log10(suplim),2)
    b = np.linspace(lim_inf, sup_lim, nn)
    x_log = 10.0**b
    fx_log = np.zeros(np.size(x_log))
    for j in range(1, np.size(x_log)-1):
        for i in range(np.size(x)-1):
            if (x_log[j] - x[i])*(x_log[j] - x[i+1]) < 0.0 :  #change of sign
                if (x_log[j] - x[i]) < (x_log[j] - x[i+1]):
                    x_log[j] = x[i]
                    fx_log[j] = f_x[i]
                else:
                    x_log[j] = x[i+1]
                    fx_log[j] = f_x[i+1]
    if suplim ==1 :    
        x_log[np.size(x_log)-1] = x[np.size(x)-1]        
        fx_log[np.size(x_log)-1] = f_x[np.size(x)-1]
    else:
        x_log[np.size(x_log)-1] = x[int(suplim/av_dt(x))]
        fx_log[np.size(x_log)-1] = f_x[int(suplim/av_dt(x))]
    if liminf == 1:
        x_log[0] = x[0]        
        fx_log[0] = f_x[0]
    else:
        x_log[0] = x[int(liminf/av_dt(x))]
        fx_log[0] = f_x[int(liminf/av_dt(x))]
    return x_log, fx_log

def linear_fit(x, y):
    """This function receives a function and performs linear fit"""
    """Input: array of dependent variable, x. And array of dependent variable, y."""
    """Output: slope, intersect, and coefficient of determination (r^2 value) """
    """An alternative is to use: slope, intercept, r_value, p_value, std_err = stats.linregress(t, defl)"""
    m,b = np.polyfit(x, y, 1)
    mean = sum(y)/np.size(y)
    SS_tot = sum((y-mean)**2)
    SS_res = sum(   (y - (m*x+b))**2     )
    r_2 = 1.0 - SS_res/SS_tot
    return m, b, r_2

def linear_fit_Nob(x,y):
    m, b, _ = linear_fit(x,y)
    params = Parameters()
    params.add('A', value = m, min=0)
    result = minimize(linear_force, params, args=(x,y), method='leastsq')
    Fdot = result.params['A'].value  
    return Fdot
    

def concatenate(t1, t2, f1, f2):
    """This function concatenates two time functions in a single one"""
    """Input: two time arrays and two time dependent functions"""
    """Output: one time array with the range including both time arrays, and two independent function corresponding to that joined time array"""
    t3 = np.zeros(np.size(t1)+np.size(t2))
    f3 = np.zeros(np.size(f1)+np.size(f2))
    for i in range(np.size(t3)):
        if i < np.size(t1):
            t3[i] = t1[i]
            f3[i] = f1[i]
        else:
            t3[i] = t2[i-np.size(t1)]
            f3[i] = f2[i-np.size(t1)]
    return t3, f3

def func_chi(params, t, chi_exp, arms=3):
    p = params.valuesdict()
    Jg = p['Jg']
    J1 = p['J1']
    tau1 = p['tau1']
    if arms > 1:
        J2 = p['J2']
        tau2 = p['tau2']
        if arms > 2:
            J3 = p['J3']
            tau3 = p['tau3']
            if arms >3:
                J4 = p['J4']
                tau4 = p['tau4']
                if arms >4:
                    J5 = p['J5']
                    tau5 = p['tau5']
                else:
                    J5 = 0.0
                    tau5 = 0.0
            else:
                J4 = 0.0
                tau4 = 0.0
                J5 = 0.0
                tau5 = 0.0
        else:
            J3 = 0.0
            tau3 = 0.0
            J4 = 0.0
            tau4 = 0.0
            J5 = 0.0
            tau5 = 0.0
    else:
        J2 = 0.0
        tau2 = 0.0
        J3 = 0.0
        tau3 = 0.0
        J4 = 0.0
        tau4 = 0.0    
        J5 = 0.0
        tau5 = 0.0           
    model = (Jg+J1+J2+J3+J4+J5)*t + J1*tau1*(exp(-t/tau1)-1.0) \
        + J2*tau2*(exp(-t/tau2)-1.0)  + J3*tau3*(exp(-t/tau3)-1.0) + J4*tau4*(exp(-t/tau4)-1.0) \
        + J5*tau5*(exp(-t/tau5)-1.0) #+ J6*tau_c[5]*(exp(-t/tau_c[5])-1.0)
    return  (model - chi_exp) /chi_exp #calculating the residual        



def func_conv(params, t, F, t_res_fit, t_exp, tip_norm_log, arms, dt):
    """This function performs a non-linear fit based on Eq. 13 and 16 in: López‐Guerra, Enrique A., Babak Eslami, and Santiago D. Solares. Journal of Polymer Science Part B: Polymer Physics (2017)."""
    p = params.valuesdict()
    Jg = p['Jg']
    J1 = p['J1']
    tau1 = p['tau1']
    U1 = J1/tau1*np.exp(-t/tau1)
    if arms > 1:
        J2 = p['J2']
        tau2 = p['tau2']
        U2 = J2/tau2*np.exp(-t/tau2)
        if arms > 2:
            J3 = p['J3']
            tau3 = p['tau3']
            U3 = J3/tau3*np.exp(-t/tau3)
            if arms >3:
                J4 = p['J4']
                tau4 = p['tau4']
                U4 = J4/tau4*np.exp(-t/tau4)
                if arms >4:
                    J5 = p['J5']
                    tau5 = p['tau5']
                    U5 = J5/tau5*np.exp(-t/tau5)
                else:
                    J5 = 0.0
                    tau5 = 0.0
                    U5 = 0.0
            else:
                J4 = 0.0
                tau4 = 0.0
                U4 = 0.0
                J5 = 0.0
                tau5 = 0.0
                U5 = 0.0
        else:
            J3 = 0.0
            tau3 = 0.0
            U3 = 0.0
            J4 = 0.0
            tau4 = 0.0
            U4 = 0.0
            J5 = 0.0
            tau5 = 0.0
            U5 = 0.0
    else:
        J2 = 0.0
        tau2 = 0.0
        U2 = 0.0
        J3 = 0.0
        tau3 = 0.0
        U3 = 0.0
        J4 = 0.0
        tau4 = 0.0    
        U4 = 0.0
        J5 = 0.0
        tau5 = 0.0    
        U5 = 0.0
    U_t = U1+U2+U3+U4+U5
    model = np.convolve(U_t, F,mode='full')*dt
    model =  model[range(np.size(F))] + Jg*F 
    model_log, _ = log_scale(model, t, t_res_fit, t_exp)
    return  (model_log - tip_norm_log) /tip_norm_log #calculating the residual
       

def chi_fit(t_simul, tip_simul, F_simul, R, t_res, t_exp, arms =3, technique=0, Jg_i = 2.0e-10, J1_i = 1.0e-9, tau1_i = 1.0e-4, J2_i = 1.0e-8, tau2_i = 1.0e-3, J3_i = 1.0e-7, tau3_i = 1.0e-2, J4_i = 1.0e-6, tau4_i = 1.0e-1, J5_i = 1.0e-5, tau5_i = 1.0e0):
    """This function performs a non-linear fit based on Eq. 14 and 15 in: López‐Guerra, Enrique A., Babak Eslami, and Santiago D. Solares. Journal of Polymer Science Part B: Polymer Physics (2017)."""
    # create a set of Parameters
    params = Parameters()
    params.add('Jg', value = Jg_i, min=0)
    params.add('J1', value = J1_i, min=0)
    params.add('tau1', value = tau1_i, min=tau1_i/10.0, max=tau1_i*10.0)
    if arms > 1:
        params.add('J2', value = J2_i, min=0)
        params.add('tau2', value = tau2_i, min=tau2_i/10.0, max=tau2_i*10.0)
        if arms >2:
            params.add('J3', value = J3_i, min=0)
            params.add('tau3', value = tau3_i, min=tau3_i/10.0, max=tau3_i*10.0)
            if arms>3:
                params.add('J4', value = J4_i, min=0)
                params.add('tau4', value = tau4_i, min=tau4_i/10.0, max=tau4_i*10.0)
                if arms>4:
                    params.add('J5', value = J5_i, min=0)
                    params.add('tau5', value = tau5_i, min=tau5_i/10.0, max=tau5_i*10.0)
    if technique ==0:
        F_log, t_log = log_scale(F_simul, t_simul, t_res, t_exp)    
        Fdot = linear_fit_Nob(t_log, F_log)    
        chi_exp = 16.0/3*np.sqrt(R)*tip_simul**1.5/Fdot    
        chi_exp_log, t_log = log_scale(chi_exp, t_simul, t_res, t_exp)
        result = minimize(func_chi, params, args=(t_log, chi_exp_log, arms), method='leastsq')
    else:
        F, t = smear(F_simul, t_simul, t_res, t_exp)  #this arrays will be passed on the convolution and have to have the time resolution of the experiment
        tip, _ = smear(tip_simul, t_simul, t_res, t_exp)
        t_res_fit = t_res*1.0        
        tip_norm = 16.0/3*np.sqrt(R)*tip**1.5
        tip_norm_log, _ = log_scale(tip_norm, t, t_res_fit, t_exp)
        dt = av_dt(t)
        result = minimize(func_conv, params, args=(t, F, t_res_fit, t_exp, tip_norm_log, arms, dt), method='leastsq')
    Jg_c = result.params['Jg'].value
    J_c= np.zeros(arms)  #N is the number of voigt units retrieved
    tau_c = np.zeros(arms)
    J_c[0] = result.params['J1'].value
    tau_c[0]= result.params['tau1'].value
    if arms > 1:
        J_c[1] =result.params['J2'].value
        tau_c[1]= result.params['tau2'].value
        if arms >2:
            J_c[2] =result.params['J3'].value
            tau_c[2]= result.params['tau3'].value
            if arms>3:
                J_c[3] =result.params['J4'].value
                tau_c[3]= result.params['tau4'].value
                if arms>4:
                    J_c[4] =result.params['J5'].value
                    tau_c[4]= result.params['tau5'].value
    return Jg_c, tau_c, J_c

    
def position(array, value):
    """finds the position of the array that contains a given value. If the value is not there it will return the nearest position"""
    """Input: array and given value"""
    """Output: nearest array position containing the value given"""
    if array[2] > array[1]:
        a = True
    else:
        a = False
    i=0
    pos=0
    flag = False
    if a:
        for i in range(np.size(array)):
            if (array[i] > value) and flag==False:
                flag = True
                if abs(array[i] - value) < abs(array[i-1]-value):
                    pos=i
                else:
                    pos= i-1
    else:
        for i in range(np.size(array)):
            if (array[i] < value) and flag==False:
                flag = True
                if abs(array[i] - value) < abs(array[i-1]-value):
                    pos=i
                else:
                    pos= i-1
    return pos
    

def plot_lims(x_array, xlim_inf, xlim_sup, array_1, array_2):
    """this function returns the inferior limit and superior limit to plot two arraysgiven certain x limits"""
    if array_1[position(x_array,xlim_inf)] < array_1[position(x_array,xlim_sup)]:  #determining if the functions to plot are growing or decreasing
        growing = True
    else:
        growing = False
    if growing:
        if array_1[position(x_array,xlim_inf)] < array_2[position(x_array,xlim_inf)]:
            c = array_1[position(x_array,xlim_inf)]
        else:
            c = array_2[position(x_array,xlim_inf)]
            
        if array_1[position(x_array,xlim_sup)] > array_2[position(x_array,xlim_sup)]:
            d = array_1[position(x_array,xlim_sup)]
        else:
            d = array_2[position(x_array,xlim_sup)]
    else:
        if array_1[position(x_array,xlim_inf)] > array_2[position(x_array,xlim_inf)]:
            d = array_1[position(x_array,xlim_inf)]
        else:
            d = array_2[position(x_array,xlim_inf)]
            
        if array_1[position(x_array,xlim_sup)] < array_2[position(x_array,xlim_sup)]:
            c = array_1[position(x_array,xlim_sup)]
        else:
            c = array_2[position(x_array,xlim_sup)]
    return c, d
    
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
