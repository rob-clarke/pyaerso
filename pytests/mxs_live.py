#!/usr/bin/env python3

import math

from numpy import mean, var

import pyaerso
from pyaerso import AffectedBody, AeroEffect, AeroBody, Body, Force, Torque

mass = 2.0
inertia = [
    [0.05, 0.0, 0.0],
    [0.0, 0.05, 0.0],
    [0.0, 0.0, 0.07]]
position = [0.0,0.0,-500.0]
velocity = [13.0,0.0,0.0]
attitude = [0.0,0.0,0.0,1.0]
rates = [0.0,0.0,0.0]

class WindModel:
    def get_wind(self,position):
        return [0,0,0]
    def step(self,dt):
        pass

import isacalc

class DensityModel:
    def get_density(self,position):
        [_,_,density,_,_] = isacalc.calculate_at_h(-position[2])
        return density

def H(x,p,k=10):
    return 1.0 / ( 1.0 + math.exp(-2*k*(x-p)) )

def S(x,l,h,k=10):
    return H(x,l,k)*(1-H(x,h,k))

def deg2rad(deg):
    return deg/180.0 * math.pi

Sw = 0.263
cw = 0.24
S_t = 0.0825
x_t = -0.585
c_t = 0.165

class Lift:
    def get_effect(self,airstate,rates,input):
        alpha = airstate[0]
        c_l = S(alpha,-0.34,0.2719) * (0.1615 + 5.2212*alpha)
        
        q = airstate[3]
        lift = q * Sw * c_l
        
        return (
            Force.body([0,0,-lift]),
            Torque.body([0,0,0])
            )

class Drag:
    def get_effect(self,airstate,rates,input):
        alpha = airstate[0]
        a_lim = deg2rad(30)
        c_d = S(alpha,-a_lim,a_lim)*(2.3814*(alpha-0.0207)**2+0.0671) + 2*(1-S(alpha,-a_lim,a_lim))
        
        q = airstate[3]
        drag = q * Sw * c_d
        
        return (
            Force.body([-drag,0,0]),
            Torque.body([0,0,0])
            )

class Moment:
    def get_effect(self,airstate,rates,input):
        alpha = airstate[0]
        a_lim = deg2rad(15)
        c_m = S(alpha,-a_lim,a_lim,12) * (-0.5462 * math.tan(1.4151*(alpha-0.0484)) + 0.053) \
            + 0.5 * (1-H(alpha,-a_lim,12)) - 0.5 * H(alpha,a_lim,12)
        
        q = airstate[3]
        V = airstate[2]
        
        alpha_q = math.atan( (-rates[1]*x_t + V*math.sin(alpha)) / (V*math.cos(alpha)) )
        m_t_aq = x_t * q * S_t * (c_lta(alpha_q) - c_lta(alpha))
        
        moment = q * Sw * cw * c_m + m_t_aq
        
        return (
            Force.body([0,0,0]),
            Torque.body([0,moment,0])
            )

def c_lta(alpha):
    return 3.5810*S(alpha,-0.1745,0.1745)*alpha \
        + 0.65 * H(alpha,0.1745) \
        - 0.65 * (1-H(alpha,-0.1745))
    

def poly_manifold(x,y,z,order,*args):    
    if order == 0:
        return args[0]
    if order == 1:
        return args[0] * x + args[1] * y + args[2] * z + args[3]
    if order == 2:
        return args[0]*x*x + args[1]*x*y + args[2]*x*z + args[3]*y*y + args[4]*y*z + args[5]*z*z \
            + args[6]*x + args[7]*y + args[8]*z \
            + args[9]
    if order == 3:
        return args[0]*x*x*x + args[1]*x*x*y + args[2]*x*x*z + args[3]*x*y*y + args[4]*x*y*z + args[5]*x*z*z + args[6]*y*y*x + args[7]*y*y*y + args[8]*y*y*z + args[9]*y*z*z + args[10]*z*z*z \
            + args[11]*x*x + args[12]*x*y + args[13]*y*y + args[14]*x*z + args[15]*y*z + args[16]*z*z \
            + args[17]*x + args[18]*y + args[19]*z \
            + args[20]

def clamp(x,u,l):
    return max(min(x,u),l)

class Combined:
    def __init__(self):
        self.time = 0

    def get_elevator_moment_coeff_delta(self,airstate,rates,input):
        coeffs = [
            -1.09419902e+00,
            -2.11376232e-01,
            2.68838823e-02,
            -2.19078045e+02,
            -1.85268630e-01,
            4.81678617e-03,
            2.19458992e+02,
            -4.51607587e-01,
            2.38686253e-03,
            -1.87506464e-04,
            1.22846953e-03,
            -4.61863302e-01,
            3.63517580e+00,
            5.27810084e-01,
            -1.51973289e-01,
            4.08504874e-03,
            -5.98825396e-02,
            2.59581816e+00,
            -1.73627947e-01,
            9.31570574e-01,
            -4.57208842e+00
        ]
        elev = input[1]
        throttle = input[2]
        aspd = clamp(airstate[2],10,22)
        return poly_manifold(elev,throttle,aspd,3,*coeffs)
    
    def get_thrust(self,airspeed,throttle):
        return 1.1576e1 * throttle**2 \
            + -4.8042e-1 * throttle * airspeed \
            + -1.1822e-2 * airspeed**2 \
            + 1.3490e1 * throttle \
            + 4.6518e-1 * airspeed \
            + -4.1090
    
    def get_effect(self,airstate,rates,input):
        alpha = airstate[0]
        c_l = S(alpha,-0.34,0.2719) * (0.1615 + 5.2212*alpha)
                
        a_lim = deg2rad(30)
        c_d = S(alpha,-a_lim,a_lim)*(2.3814*(alpha-0.0207)**2+0.0671) + 2*(1-S(alpha,-a_lim,a_lim))
        
        a_lim = deg2rad(15)
        c_m = S(alpha,-a_lim,a_lim,12) * (-0.5462 * math.tan(1.4151*(alpha-0.0484)) + 0.053) \
            + 0.5 * (1-H(alpha,-a_lim,12)) - 0.5 * H(alpha,a_lim,12)
        
        c_m += x_t * c_lta(alpha) / c_t # * (1-S(alpha,-a_lim,a_lim,12))
        
        q = airstate[3]
        V = airstate[2]
        
        alpha_q = math.atan2( (-rates[1]*x_t + V*math.sin(alpha)), (V*math.cos(alpha)) )
        m_t_aq = x_t * q * S_t * (c_lta(alpha_q) - c_lta(alpha))
        
        dc_m_elev = self.get_elevator_moment_coeff_delta(airstate,rates,input)
        
        lift = q * Sw * c_l
        drag = q * Sw * c_d
        moment = q * Sw * cw * (c_m + dc_m_elev) + m_t_aq
        
        thrust = self.get_thrust(V,input[2])
        # print(f"t: {self.time:.3f} T: {thrust:.3f}, V: {V:.3f}, Alpha: {alpha:.3f}, q: {q:.3f},  c_m: {c_m:.3f}, alpha_q: {alpha_q:.3f}, c_m_eff: {c_m+dc_m_elev:.3f}, moment: {moment:.3f}, dc_m_elev: {dc_m_elev:.3f}, mtaq: {m_t_aq:.3f}")
        self.time += 0.01
        x = thrust - drag * math.cos(alpha) + lift * math.sin(alpha)
        z = -lift * math.cos(alpha) - drag * math.sin(alpha)
        
        return (
            Force.body([x,0.0,z]),
            Torque.body([0,moment,0])
            )


def sensible(array,dp=4):
    elements = ", ".join([f"{v:.{dp}f}" for v in array])
    return f"[{elements}]"

body = Body(mass,inertia,position,velocity,attitude,rates)

# aerobody = AeroBody(body,WindModel(),("StandardDensity",[]))
# aerobody = AeroBody(body)
aerobody = AeroBody(body,None,DensityModel())

# vehicle = AffectedBody(aerobody,[Lift(),Drag(),Moment()])
vehicle = AffectedBody(aerobody,[Combined()])

# print(vehicle.airstate)

scale = 2
deltaT = 0.01/scale

print(sensible(vehicle.statevector))

import sys
import time

outfile = None
if len(sys.argv) > 1:
    outfile = open(sys.argv[1],"w")
    outfile.write("time,x,y,z,u,v,w,qx,qy,qz,qw,p,q,r,alpha\n")

samples = 1
sample_times = []

count = 0
simtime = 0

body = Body(mass,inertia,position,velocity,attitude,rates)
aerobody = AeroBody(body,None,("StandardDensity",[]))
vehicle = AffectedBody(aerobody,[Combined()])

frameInterval = 20

times = []
airspeeds = []
qs = []
xs = []
zs = []

steps_per_frame = (frameInterval/1000)/deltaT

def render(frame, *args):
    for i in range(round(steps_per_frame)):
        global count
        global simtime
        elevator = (0.0/180.0)*math.pi if count < 500*scale else ( (5.0/180.0)*math.pi if count < 900*scale else (0.0/180.0)*math.pi )
        # throttle = 0.7687
        throttle = 0.78
        vehicle.step(deltaT,[0,elevator,throttle])
        count += 1
        simtime += deltaT
    times.append(simtime)
    airspeeds.append(vehicle.airstate[2])
    xs.append(vehicle.statevector[0])
    zs.append(-vehicle.statevector[2])
    qs.append(math.degrees(vehicle.statevector[11]))
    
    [position_trace, airspeed_trace, q_trace] = args
    
    position_trace.set_xdata(xs)
    position_trace.set_ydata(zs)
    
    airspeed_trace.set_xdata(times)
    airspeed_trace.set_ydata(airspeeds)

    q_trace.set_xdata(times)
    q_trace.set_ydata(qs)

    return [position_trace, airspeed_trace, q_trace]

import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, axs = plt.subplots(2,2)
axs[0,0].set_xlim(0,200)
axs[0,0].set_ylim(200,700)
axs[0,0].axis('equal')
axs[0,0].grid('both')

axs[0,1].set_xlim(0,60)
axs[0,1].set_ylim(0,40)
axs[0,1].grid('both')

axs[1,0].set_xlim(0,60)
axs[1,0].set_ylim(-90,90)
axs[1,0].grid('both')

position_trace, = axs[0,0].plot(xs,zs)
airspeed_trace, = axs[0,1].plot(times,airspeeds)
q_trace, = axs[1,0].plot(times,qs)
ani = animation.FuncAnimation(fig, render, interval=frameInterval, fargs=[position_trace, airspeed_trace, q_trace], blit=True)

plt.show()
