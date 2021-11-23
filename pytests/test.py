#!/usr/bin/env python3

import pyaerso
from pyaerso import AeroBody, Body, Force

mass = 1.0
inertia = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]]
position = [0.0,0.0,0.0]
velocity = [0.0,0.0,0.0]
attitude = [0.0,0.0,0.0,1.0]
rates = [0.0,0.0,0.0]

class WindModel:
    def get_wind(self,position):
        return [0,0,0]
    def step(self,dt):
        pass

wm = WindModel()

body = Body(mass,inertia,position,velocity,attitude,rates)
vehicle = AeroBody.new(body,wm,("StandardDensity",[]))

vehicle.airstate

thrust = Force.body([4.0,0.0,0.0])

def get_thrust(power):
    return -0.0000830488*power**2 + 0.0704307060*power + 0.5996810096

print(vehicle.attitude)

deltaT = 0.01
time = 0
while vehicle.position[2] < 2000:
    vehicle.step([thrust],[],deltaT)
    time += deltaT

print(time)
print(vehicle.statevector)
