#!/usr/bin/env python3

import pyaerso
from pyaerso import AffectedBody, AeroEffect, AeroBody, Body, Force, Torque

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

class Lift:
    def get_effect(self,airstate,rates,input):
        print(f"airstate: {airstate}")
        print(f"rates: {rates}")
        print(f"input: {input}")
        return (
            Force.body([0,0,0]),
            Torque.body([0,0,0])
            )

class Thrust:
    def get_thrust(self,power):
        return -0.0000830488*power**2 + 0.0704307060*power + 0.5996810096
    
    def get_effect(self,airstate,rates,input):
        thrust = self.get_thrust(input[2])
        return (
            Force.body([thrust,0,0]),
            Torque.body([0,0,0])
            )


body = Body(mass,inertia,position,velocity,attitude,rates)
aerobody = AeroBody(body,wm,("StandardDensity",[]))
vehicle = AffectedBody(aerobody,[Lift(),Thrust()])

vehicle.airstate

print(vehicle.attitude)

deltaT = 0.01
time = 0
while vehicle.position[2] < 2000:
    vehicle.step(deltaT,[0,0,100])
    time += deltaT

print(time)
print(vehicle.statevector)
