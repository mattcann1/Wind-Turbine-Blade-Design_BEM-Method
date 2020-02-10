# -----------------------------------------------------------#
# (C) 2020 Matthew Cann
# Released under MIT Public License (MIT)
# email mcann@uwaterloo.ca
# -----------------------------------------------------------

'''Performs BEM method of wind turbine design using the operating conditions
of a Vestas V90-2MW wind turbine.'''

#IMPORTS.......................................................................
import math
import numpy as np

from pandas import DataFrame
from scipy.stats import exponweib 
from scipy.stats import rayleigh
import scipy
import scipy.stats
import pandas as pd
import matplotlib.pyplot as plt

#CONSTANTS.....................................................................
UU = 13.5 #m/s          #Rated Wind Speed of Vestas V90-2MW Wind Turbine
RR = 45 #m              #Radius of Vestas V90-2MW Wind Turbine swept area
rho = 1.225 #kg/m3      #Density of Air
Power = 2 *10**6 #W     #Rated Power of Vestas V90-2MW Wind Turbine
Rot_speed = 14.9#rpm    #Rotation Speed of Vestas V90-2MW Wind Turbine
mu = 0.00001802         #Air Mu 
E = 20 # number of blade elements
hub = 0.1 # perctage of the blade that the hub takes up
N = 3 #Blades

#DEFINTIONS...................................................................

def linear_blade(blade_phi, blade_chord,  r_R):
    #Linearize
    Lin_phi = [blade_phi[-2], blade_phi[-1]]
    Lin_chord = [blade_chord[-2], blade_chord[-1]]
    
    z_phi = np.polyfit([r_R[-2], r_R[-1]], Lin_phi, 1)
    z_chord = np.polyfit([r_R[-2], r_R[-1]], Lin_chord, 1)
    
    lin_phi_eq = str(z_phi[0]) + "*x + " + str(z_phi[1]) # x is r/R
    lin_chord_eq = str(z_chord[0]) + "*x + " + str(z_chord[1])
    
    #print(lin_phi_eq, lin_chord_eq)
    
    lin_phi = [eval(lin_phi_eq) for x in r_R]
    lin_chord = [eval(lin_chord_eq) for x in r_R]
    
    return (lin_phi, lin_chord)

def rotor_performance(Radius,r_R,l_lamda,N,Cl_max, dr, l_chord):
    percent = 100
    count = 0
    a = 1/3
    a_p = 0    
    while percent > 0.1:
        l_radius = r_R   
        tan_phi = (1-a)/(l_lamda*(1+a_p))
        l_phi = math.degrees(math.atan(tan_phi))
        l_coef_lift = Cl_max
        sigma_r = (N*l_chord)/(2*math.pi*l_radius*Radius)
        coefa = (sigma_r*l_coef_lift*math.cos(math.radians(l_phi)))/(4*(math.sin(math.radians(l_phi)))**2)
        #coefa_p = (sigma_r*l_coef_lift)/(4*math.cos(math.radians(l_phi)))
        coefa_p = (4*math.cos(math.radians(l_phi))-sigma_r*l_coef_lift)/(sigma_r*l_coef_lift)
        u_a = coefa/(1+coefa)
        #u_a_p = coefa_p/(1-coefa_p)
        u_a_p = 1/coefa_p
        percent = abs((u_a-a)/a*100)
        a = u_a
        a_p = u_a_p
        #print(count, l_phi, a, a_p)
        count+=1
    return a, a_p, l_phi

#Initial Calculations and Unit Transformations.................................
#Converts Rotatation into rad/sec
Omega = Rot_speed *2*math.pi/60 #rad/sec

#Estimate Coefficent of Performance (Cp)
eff = 0.5
Cp_est = Power/(eff*0.5*rho*math.pi*RR**2*UU**3)
print(Cp_est)

#Tip speed ratio of rated wind turbine. 
Lamda = (Omega)*RR/UU
print(Lamda)


#%%
#Airfoil Information...........................................................

#Reynolds Number 1 million
#Airfoil: EPPLER 58 AIRFOIL
Cl_1 = 1.1054 # cl where Cl/Cd is max @1.75
alpha_1 = 1.75 #degrees
Cd_1 = 0.00469

#Reynolds Number 2 million
#Airfoil: DAVIS BASIC B-24 WING AIRFOIL
alpha_2 = 7.500
Cl_2 = 1.0440
Cd_2 = 0.00759

#Reynolds Number 5 million
#Airfoil: DAVIS BASIC B-24 WING AIRFOIL
alpha_3 = 7
Cl_3 = 0.9377
Cd_3 = 0.00595

#Element Analysis..............................................................

#Optimal Parameters
a = 1/3
a_p = 0


blade_lamda = []
blade_phi = []
blade_chord = []
r_R = []

for e in range(1, E+1):
    
    dr = 1/E
    local_radius = e*dr-dr/2
    r_R.append(local_radius)
    
    local_lamda = Lamda*local_radius
    local_phi = 2/3*math.degrees(math.atan(1/local_lamda))
    if e < 7:
        local_coef_lift = Cl_1
        local_chord = (8*math.pi*local_radius*RR)/(N*local_coef_lift)*(1-math.cos(math.radians(local_phi)))
        blade_lamda.append(local_lamda)
        blade_phi.append(local_phi)
        blade_chord.append(local_chord)
        
    if 7 <= e < 14:
        local_coef_lift = Cl_2
        local_chord = (8*math.pi*local_radius*RR)/(N*local_coef_lift)*(1-math.cos(math.radians(local_phi)))
        blade_lamda.append(local_lamda)
        blade_phi.append(local_phi)
        blade_chord.append(local_chord)
      

    if e >= 14:
        local_coef_lift = Cl_3
        local_chord = (8*math.pi*local_radius*RR)/(N*local_coef_lift)*(1-math.cos(math.radians(local_phi)))
        blade_lamda.append(local_lamda)
        blade_phi.append(local_phi)
        blade_chord.append(local_chord)
      
df = DataFrame({'r/R': r_R, '\u03BB r': blade_lamda, "\u03C6 r [deg]": blade_phi, 'Chord [m]': blade_chord})

#%%
#Overal Blade shape=====================================================

lin_phi, lin_chord = linear_blade(blade_phi, blade_chord,r_R)

df['linear \u03C6 r [deg]'] = lin_phi
df['linear chord [m]'] = lin_chord 

plt.figure(1)
plt.ylabel('Phi [deg]')
plt.xlabel('r/R')
plt.xlim(0,1)
plt.plot(r_R[2:],blade_phi[2:], label = 'Optimum Angle')
plt.plot(r_R[2:],lin_phi[2:],'--', label = 'Linearized Angle')
plt.savefig('Lin_angle.png',bbox_inches="tight")
plt.legend()
plt.show()

plt.figure(2)
plt.ylabel('Chord [m]')
plt.xlabel('r/R')
plt.xlim(0,1)

plt.plot(r_R[2:],blade_chord[2:], label = 'Optimum Chord')
plt.plot(r_R[2:],lin_chord[2:],'--', label = 'Linearized Chord')
plt.legend()
plt.savefig('Lin_chord.png',bbox_inches="tight")
plt.show()

#%%
#Rotor Performance =====================================================


a_list = []
a_p_list = []

rot_vel_l = []
wind_vel_l = []
rel_vel_l = []
Power_nodrag_list = []
Power_drag_list = []
coef_perf_l = []
coef_perf_nodrag_l = []


force_lift_list = []
force_drag_list = []
force_norm_list = []
force_tang_list = []
force_torque_list = []
power_initial = []


phi_l = []

no_drag_power_l = []
        
power_l = []
FF_l = []
TL_coef_perf_l = []
TL_power_l = []

for e in range(1, E+1):
    dr = 1/E
    l_radius = e*dr-dr/2
    l_lamda = Lamda*l_radius
    l_chord = lin_chord[e-1]
    
    if e < 7:
        Cl_max = Cl_1
        Cd = Cd_1       
        
        a, a_p, phi = rotor_performance(RR, l_radius, l_lamda, N, Cl_max, dr, l_chord)
        a_list.append(a)
        a_p_list.append(a_p)
        #print(a,a_p, phi)
        phi_l.append(phi)
        rot_vel = Omega*l_radius*RR*(1+a_p)
        wind_vel = UU*(1-a)
        rel_vel = math.sqrt(rot_vel**2+wind_vel**2) #W 
        #print(UU)
        #print(rot_vel,wind_vel, rel_vel)
        rot_vel_l.append(rot_vel) 
        wind_vel_l.append(wind_vel)
        rel_vel_l.append(rel_vel) #W
        
        Reynolds_num = rho*rel_vel*l_chord/mu
        #print(Reynolds_num)
        
        
        l_phi = lin_phi[e-1]
        #l_phi = phi
        
        l_solidity = N*l_chord/(2*math.pi*l_radius*RR)
        
        l_coef = Cl_max
        #print(local_coef_lift,l_coef )
        d_coef = Cd
        #Forces =============================================================================
        #Lift Force
        force_lift = l_coef *0.5*rho*rel_vel**2*l_chord*dr
        force_lift_list.append(force_lift)
        #Drag Force
        force_drag = d_coef*0.5*rho*rel_vel**2*l_chord*dr
        force_drag_list.append(force_drag)
        #Tangential Force
        force_tang = force_lift*math.sin(math.radians(l_phi)) - force_drag*math.cos(math.radians(l_phi))
        #Normal or axial force
        
        force_norm = force_lift*math.cos(math.radians(l_phi)) + force_drag*math.sin(math.radians(l_phi))
        #print(force_tang, force_norm)
        force_tang_list.append(force_tang)
        force_norm_list.append(force_norm)
        #Torque
        force_torque = N*force_tang*l_radius*RR
        Power_initial = force_torque*Omega
        force_torque_list.append(force_torque)
        power_initial.append(Power_initial)
        
        #IMPORTAN*************************************************************
        coef_perf = (math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*(1-(force_drag/force_lift)*(1/(math.tan(math.radians(l_phi)))))*l_lamda**2
        coef_perf_l.append(coef_perf)
        #print(coef_perf)
        coef_perf_nodrag = (math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*l_lamda**2
        #print(coef_perf_nodrag)
        coef_perf_nodrag_l.append(coef_perf_nodrag)
        
        power = eff*coef_perf*0.5*rho*math.pi*RR**2*UU**3
        power_l.append(power)
        no_drag_power = eff*coef_perf_nodrag*0.5*rho*math.pi*RR**2*UU**3
        no_drag_power_l.append(no_drag_power)
            
    if 7 <= e < 14:
        Cl_max = Cl_2
        Cd = Cd_2
        
        a, a_p, phi = rotor_performance(RR, l_radius, l_lamda, N, Cl_max, dr, l_chord)
        a_list.append(a)
        a_p_list.append(a_p)
        #print(a,a_p, phi)
        phi_l.append(phi)
        rot_vel = Omega*l_radius*RR*(1+a_p)
        wind_vel = UU*(1-a)
        rel_vel = math.sqrt(rot_vel**2+wind_vel**2) #W 
        #print(UU)
        #print(rot_vel,wind_vel, rel_vel)
        rot_vel_l.append(rot_vel) 
        wind_vel_l.append(wind_vel)
        rel_vel_l.append(rel_vel) #W
        
        Reynolds_num = rho*rel_vel*l_chord/mu
        #print(Reynolds_num)
        
        
        l_phi = lin_phi[e-1]
        #l_phi = phi
        
        l_solidity = N*l_chord/(2*math.pi*l_radius*RR)
        #local_coef_lift = (4*math.sin(math.radians(l_phi)))*(math.cos(math.radians(l_phi))-l_lamda*(math.sin(math.radians(l_phi))))/(l_solidity*(math.sin(math.radians(l_phi))+l_lamda*(math.cos(math.radians(l_phi)))))
        
        l_coef = Cl_max
        #print(local_coef_lift,l_coef )
        d_coef = Cd
        #Forces =============================================================================
        #Lift Force
        force_lift = l_coef *0.5*rho*rel_vel**2*l_chord*dr
        force_lift_list.append(force_lift)
        #Drag Force
        force_drag = d_coef*0.5*rho*rel_vel**2*l_chord*dr
        force_drag_list.append(force_drag)
        #Tangential Force
        force_tang = force_lift*math.sin(math.radians(l_phi)) - force_drag*math.cos(math.radians(l_phi))
        #Normal or axial force
        force_norm = force_lift*math.cos(math.radians(l_phi)) + force_drag*math.sin(math.radians(l_phi))
        #print(force_tang, force_norm)
        force_tang_list.append(force_tang)
        force_norm_list.append(force_norm)
        #Torque
        force_torque = N*force_tang*l_radius*RR
        Power_initial = force_torque*Omega
        force_torque_list.append(force_torque)
        power_initial.append(Power_initial)
        
        
        #Final Calcs...........................................................
        coef_perf = (math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*(1-(force_drag/force_lift)*(1/(math.tan(math.radians(l_phi)))))*l_lamda**2
        coef_perf_l.append(coef_perf)
        #print(coef_perf)
        coef_perf_nodrag = (math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*l_lamda**2
        #print(coef_perf_nodrag)
        coef_perf_nodrag_l.append(coef_perf_nodrag)
        
        power = eff*coef_perf*0.5*rho*math.pi*RR**2*UU**3
        power_l.append(power)
        no_drag_power = eff*coef_perf_nodrag*0.5*rho*math.pi*RR**2*UU**3
        no_drag_power_l.append(no_drag_power)
    if e >= 14:
        Cl_max = Cl_3
        Cd = Cd_3
        
        a, a_p, phi = rotor_performance(RR, l_radius, l_lamda, N, Cl_max, dr, l_chord)
        a_list.append(a)
        a_p_list.append(a_p)
        #print(a,a_p, phi)
        phi_l.append(phi)
        rot_vel = Omega*l_radius*RR*(1+a_p)
        wind_vel = UU*(1-a)
        rel_vel = math.sqrt(rot_vel**2+wind_vel**2) #W 
        #print(UU)
        #print(rot_vel,wind_vel, rel_vel)
        rot_vel_l.append(rot_vel) 
        wind_vel_l.append(wind_vel)
        rel_vel_l.append(rel_vel) #W
        
        Reynolds_num = rho*rel_vel*l_chord/mu
        #print(Reynolds_num)
        
        l_phi = lin_phi[e-1]
        #l_phi = phi
        
        l_solidity = N*l_chord/(2*math.pi*l_radius*RR)
        
        l_coef = Cl_max
        d_coef = Cd
        #Forces ===============================================================
        #Lift Force
        force_lift = l_coef *0.5*rho*rel_vel**2*l_chord*dr
        force_lift_list.append(force_lift)
        #Drag Force
        force_drag = d_coef*0.5*rho*rel_vel**2*l_chord*dr
        force_drag_list.append(force_drag)
        #Tangential Force
        force_tang = force_lift*math.sin(math.radians(l_phi)) - force_drag*math.cos(math.radians(l_phi))
        #Normal or axial force
        force_norm = force_lift*math.cos(math.radians(l_phi)) + force_drag*math.sin(math.radians(l_phi))
        #print(force_tang, force_norm)
        force_tang_list.append(force_tang)
        force_norm_list.append(force_norm)
        #Torque
        force_torque = N*force_tang*l_radius*RR
        
        Power_initial = force_torque*Omega
        force_torque_list.append(force_torque)
        power_initial.append(Power_initial)
        
        #IMPORTAN*************************************************************
        coef_perf = (math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*(1-(force_drag/force_lift)*(1/(math.tan(math.radians(l_phi)))))*l_lamda**2
        coef_perf_l.append(coef_perf)
        #print(coef_perf)
        coef_perf_nodrag = (math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*l_lamda**2
        #print(coef_perf_nodrag)
        coef_perf_nodrag_l.append(coef_perf_nodrag)
        
        no_drag_power = eff*coef_perf_nodrag*0.5*rho*math.pi*RR**2*UU**3
        no_drag_power_l.append(no_drag_power)
        
        power = eff*coef_perf*0.5*rho*math.pi*RR**2*UU**3
        power_l.append(power)
        
        
#Including Tiploses............................................................

def rotor_performance_TL(Radius,r_R,l_lamda,N,Cl_max, dr, l_chord,FF):
    percent = 100
    count = 0
    a = 1/3
    a_p = 0    
    while percent > 0.1:
        l_radius = r_R   
        tan_phi = (1-a)/(l_lamda*(1+a_p))
        l_phi = math.degrees(math.atan(tan_phi))
        l_coef_lift = Cl_max
        
       
   
        sigma_r = (N*l_chord)/(2*math.pi*l_radius*Radius)
        coefa = (sigma_r*l_coef_lift*math.cos(math.radians(l_phi)))/(4*FF*(math.sin(math.radians(l_phi)))**2)
        #coefa_p = (sigma_r*l_coef_lift)/(4*FF*math.cos(math.radians(l_phi)))
        coefa_p = (4*math.cos(math.radians(l_phi))-sigma_r*l_coef_lift)/(sigma_r*l_coef_lift)
        u_a = coefa/(1+coefa)
        u_a_p = 1/coefa_p
        percent = abs((u_a-a)/a*100)
        a = u_a
        a_p = u_a_p
        #print(count, l_phi, a, a_p)
        count+=1
    return a, a_p, l_phi

    
for e in range(1, E+1):
    dr = 1/E
    l_radius = e*dr-dr/2
    l_lamda = Lamda*l_radius
    l_chord = lin_chord[e-1]
    l_phi = lin_phi[e-1]
        
    ff = (N/2)*(RR-l_radius*RR)/(l_radius*RR*(math.sin(math.radians(l_phi))))
    FF = (2/math.pi)*math.acos(math.exp(-ff))
    FF_l.append(FF)

    if e < 7:
        l_coef = Cl_1
        d_coef = Cd_1
   
        a, a_p, phi = rotor_performance_TL(RR, l_radius, l_lamda, N, l_coef, dr, l_chord,FF)
        #print(a, a_p,FF)
        rot_vel = Omega*l_radius*RR*(1+a_p)
        wind_vel = UU*(1-a)
        rel_vel = math.sqrt(rot_vel**2+wind_vel**2) #W 
        #print(UU)

        #Lift Force
        force_lift = l_coef *0.5*rho*rel_vel**2*l_chord*dr
        #force_lift_list.append(force_lift)
        #Drag Force
        force_drag = d_coef*0.5*rho*rel_vel**2*l_chord*dr
        #force_drag_list.append(force_drag)
        
        TL_coef_perf = FF*(math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*(1-(force_drag/force_lift)*(1/(math.tan(math.radians(l_phi)))))*l_lamda**2
        TL_coef_perf_l.append(TL_coef_perf)
        
        TL_power = eff*TL_coef_perf*0.5*rho*math.pi*RR**2*UU**3
        TL_power_l.append(TL_power)

        
    if 7 <= e < 14:
        l_coef = Cl_2
        d_coef = Cd_2
   
        a, a_p, phi = rotor_performance_TL(RR, l_radius, l_lamda, N, l_coef, dr, l_chord,FF)
        #a, a_p, phi = rotor_performance(RR, l_radius, l_lamda, N, Cl_max, dr, l_chord)
        #print(a, a_p,FF)
        rot_vel = Omega*l_radius*RR*(1+a_p)
        wind_vel = UU*(1-a)
        rel_vel = math.sqrt(rot_vel**2+wind_vel**2) #W 
        #print(UU)
        #print(rot_vel,wind_vel, rel_vel)
        #rot_vel_l.append(rot_vel) 
        #wind_vel_l.append(wind_vel)
        #rel_vel_l.append(rel_vel) #W
        #Lift Force
        force_lift = l_coef *0.5*rho*rel_vel**2*l_chord*dr
        #force_lift_list.append(force_lift)
        #Drag Force
        force_drag = d_coef*0.5*rho*rel_vel**2*l_chord*dr
        #force_drag_list.append(force_drag)
        
        TL_coef_perf = FF*(math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*(1-(force_drag/force_lift)*(1/(math.tan(math.radians(l_phi)))))*l_lamda**2
        TL_coef_perf_l.append(TL_coef_perf)
        
        TL_power = eff*TL_coef_perf*0.5*rho*math.pi*RR**2*UU**3
        TL_power_l.append(TL_power)

    if e >= 14:
        l_coef = Cl_3
        d_coef = Cd_3
   
        a, a_p, phi = rotor_performance_TL(RR, l_radius, l_lamda, N, l_coef, dr, l_chord,FF)
        #a, a_p, phi = rotor_performance(RR, l_radius, l_lamda, N, Cl_max, dr, l_chord)
        #print(a, a_p,FF)
        rot_vel = Omega*l_radius*RR*(1+a_p)
        wind_vel = UU*(1-a)
        rel_vel = math.sqrt(rot_vel**2+wind_vel**2) #W 
        #print(UU)
        #print(rot_vel,wind_vel, rel_vel)
        #rot_vel_l.append(rot_vel) 
        #wind_vel_l.append(wind_vel)
        #rel_vel_l.append(rel_vel) #W
        #Lift Force
        force_lift = l_coef *0.5*rho*rel_vel**2*l_chord*dr
        #force_lift_list.append(force_lift)
        #Drag Force
        force_drag = d_coef*0.5*rho*rel_vel**2*l_chord*dr
        #force_drag_list.append(force_drag)
        
        TL_coef_perf = FF*(math.sin(math.radians(l_phi)))**2*(math.cos(math.radians(l_phi))-l_lamda*math.sin(math.radians(l_phi)))*(math.sin(math.radians(l_phi))+l_lamda*math.cos(math.radians(l_phi)))*(1-(force_drag/force_lift)*(1/(math.tan(math.radians(l_phi)))))*l_lamda**2
        TL_coef_perf_l.append(TL_coef_perf)
        
        TL_power = eff*TL_coef_perf*0.5*rho*math.pi*RR**2*UU**3
        TL_power_l.append(TL_power)
        
        
df['Rotational Velocity']= rot_vel_l
df['Wind Velocity']= wind_vel_l
df['Realtive Velocity W'] = rel_vel_l
df['Lift Force [N]'] = force_lift_list
df['Drag Force [N]'] = force_drag_list
df['Normal Force [N]'] = force_norm_list
df['Tangential Force [N]'] = force_tang_list
df['Force Torque']= force_torque_list 
df['Power']=power_initial
        
        
df['a'] = a_list
df["a'"] = a_p_list

df['Cp'] = coef_perf_l
df['Cp nodrag']= coef_perf_nodrag_l


df['F'] = FF_l
df['No Drag Power [W]']=no_drag_power_l
df['Drag Power [W]'] = power_l
df['Drag with TL Power [W]'] = TL_power_l

total_power_drag = sum(Power_drag_list)
total_power_nodrag = sum(Power_nodrag_list)

AA = 6362 #Swept area [m^2]

#print(power_l)
plt.figure(3)
plt.plot(r_R, power_l, '.',linestyle='-', linewidth = 1, markersize = 10)
plt.xlabel('Non-dimensional radius [r/R]')
plt.ylabel('Power [W]')
plt.savefig('Power_blade.png',bbox_inches="tight")

plt.show()


plt.figure(4)
plt.plot(r_R[2:], power_l[2:], '.',linestyle='-', linewidth = 1, markersize = 10, label = 'Drag with no tip loss')
plt.plot(r_R[2:],TL_power_l[2:],'.',linestyle='-', linewidth = 1, markersize = 10, label = 'Drag with tip loss')
plt.legend()
plt.xlabel('Non-dimensional radius [r/R]')
plt.xlim(0,1)
plt.ylabel('Power [W]')
plt.savefig('Power_comp.png',bbox_inches="tight")
plt.show()


TL_coef_performance = sum(TL_coef_perf_l)*8/(Lamda*E)
coef_performance_drag = sum(coef_perf_l)*8/(Lamda*E)
coef_performance_nodrag = sum(coef_perf_nodrag_l)*8/(Lamda*E)

Power_drag = eff*coef_performance_drag*0.5*rho*math.pi*RR**2*UU**3/1000000
Power_nodrag = eff*coef_performance_nodrag*0.5*rho*math.pi*RR**2*UU**3/1000000
Power_TL = eff*TL_coef_performance*0.5*rho*math.pi*RR**2*UU**3/1000000

print('Cp %0.3f' %coef_performance_drag)
print('Cp - no drag %0.3f' %coef_performance_nodrag)
print('Tip loss Cp %0.3f' %TL_coef_performance)
print('Power no drag %4.3f' %Power_nodrag,'[MW]')
print('Power blade drag %4.3f' %Power_drag, '[MW]')
print('Power TL %4.3f' %Power_TL,'[MW]')



#print(df)
#df.to_excel('test5.xlsx', sheet_name='sheet1', index=False)

#%%



