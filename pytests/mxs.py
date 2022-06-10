#!/usr/bin/env python3

import math

from numpy import mean, var

import pyaerso
from pyaerso import AffectedBody, AeroEffect, AeroBody, Body, Force, Torque

mass = 2.0
inertia = [
    [0.1, 0.0, 0.0],
    [0.0, 0.1, 0.0],
    [0.0, 0.0, 0.1]]
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
    
class Combined:
    def get_effect(self,airstate,rates,input):
        alpha = airstate[0]
        c_l = S(alpha,-0.34,0.2719) * (0.1615 + 5.2212*alpha)
                
        a_lim = deg2rad(30)
        c_d = S(alpha,-a_lim,a_lim)*(2.3814*(alpha-0.0207)**2+0.0671) + 2*(1-S(alpha,-a_lim,a_lim))
        
        a_lim = deg2rad(15)
        c_m = S(alpha,-a_lim,a_lim,12) * (-0.5462 * math.tan(1.4151*(alpha-0.0484)) + 0.053) \
            + 0.5 * (1-H(alpha,-a_lim,12)) - 0.5 * H(alpha,a_lim,12)
        
        q = airstate[3]
        V = airstate[2]
        
        alpha_q = math.atan( (-rates[1]*x_t + V*math.sin(alpha)) / (V*math.cos(alpha)) )
        m_t_aq = x_t * q * S_t * (c_lta(alpha_q) - c_lta(alpha))
        
        lift = q * Sw * c_l
        drag = q * Sw * c_d
        moment = q * Sw * cw * c_m + m_t_aq
        
        # print(f"Alpha: {alpha:.3f}, q: {q:.3f},  c_m: {c_m:.3f}, alpha_q: {alpha_q:.3f}, moment: {moment:.3f}")
        
        return (
            Force.body([-drag,0,-lift]),
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

scale = 1
deltaT = 0.01/scale

print(sensible(vehicle.statevector))

import sys
import time

outfile = None
if len(sys.argv) > 1:
    outfile = open(sys.argv[1],"w")
    outfile.write("time,x,y,z,u,v,w,qx,qy,qz,qw,p,q,r\n")

samples = 5
sample_times = []

for i in range(samples):
    count = 0
    simtime = 0
    
    body = Body(mass,inertia,position,velocity,attitude,rates)
    # aerobody = AeroBody(body,None,DensityModel())
    aerobody = AeroBody(body,None,("StandardDensity",[]))
    vehicle = AffectedBody(aerobody,[Combined()])
    
    start = time.process_time()
    while count < 6000*scale:
        vehicle.step(deltaT,[0,0,0])
        count += 1
        simtime += deltaT
        vehicle.statevector
        if outfile:
            outfile.write(f"{simtime},"+sensible(vehicle.statevector,10)[1:-1] + "\n")
    end = time.process_time()
    sample_times.append(end-start)
    print(sensible(vehicle.statevector))

print(f"Mean: {mean(sample_times)}\nVar: {var(sample_times)}")
if outfile is not None:
    outfile.close()
