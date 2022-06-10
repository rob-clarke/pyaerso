# `pyaerso`

`pyaerso` is a Python wrapper around `aerso`, a Rust library for kinematics with helpers for aerodynamic forces.

## Example

> The full code for the example can be found in [`pytests/test.py`](./blob/main/pytests/test.py).

The first step is to create a `Body` which represents the base kinematic properties:

```py
from pyaerso import AffectedBody, AeroBody, Body, Force, Torque

mass = 1.0
inertia = [
    [0.1, 0.0, 0.0],
    [0.0, 0.1, 0.0],
    [0.0, 0.0, 0.1]]
position = [0.0,0.0,0.0]
velocity = [0.0,0.0,0.0]
attitude = [0.0,0.0,0.0,1.0]
rates = [0.0,0.0,0.0]

body = Body(mass,inertia,position,velocity,attitude,rates)
```

This is then wrapped in an `AeroBody`, where wind and density models can be attached:

```py
class WindModel:
    def get_wind(self,position):
        return [0,0,0]
    def step(self,dt):
        pass

aerobody = AeroBody(body,WindModel(),("StandardDensity",[]))
```

Next any `AeroEffect`s can be defined:

```py
class Lift:
    def get_effect(self,airstate,rates,input):
        alpha = (airstate[0]/math.pi)*180.0
        q = airstate[3]
        c_l = -(1.0/6500.0)*alpha**3 + 0.1*alpha if abs(alpha) < 20.0 else math.copysign(0.5,alpha)
        lift = q * 0.3 * c_l
        return (
            Force.body([0,0,-lift]),
            Torque.body([0,0,0])
            )

class Drag:
    def get_effect(self,airstate,rates,input):
        alpha = (airstate[0]/math.pi)*180.0
        q = airstate[3]
        c_d = (1/1300.0)*(alpha-2)**2 + 0.07 if abs(alpha) < 20.0 else 0.3
        drag = q * 0.3 * c_d
        return (
            Force.body([-drag,0,0]),
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
```

Finally, the effects are attached via the `AffectedBody`:

```py
vehicle = AffectedBody(aerobody,[Lift(),Thrust(),Drag()])
```

Then propogating the system is just calling `step`:

```py
deltaT = 0.01
time = 0
time_limit = 5*2000
while time < time_limit:
    vehicle.step(deltaT,[0,0,200])
    time += deltaT
```

## `pyaerso_py`

This is a "pure" Python implementation of the module for comparison and educational purposes.

There is likely some further testing to be done.

The interface is *(hopefully)* identical to the Rust-backed version.
