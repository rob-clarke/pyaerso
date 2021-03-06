from abc import ABC, abstractmethod
from typing import List, Tuple
import numpy as np

from .types import Force, Torque

from .Aero import AeroBody, AirState

class AeroEffect(ABC):
    """Trait for aerodynamic effect"""
    
    @abstractmethod
    def get_effect(self, airstate: AirState, rates: List[float], inputstate: List[float]) -> Tuple[Force,Torque]:
        """
        Return `Force` and `Torque` generated by the effect.
        
        When this effect is given to an [AffectedBody], `get_effect` will be called at each timestep
        to get the force and torque induced by the effect. The returned Force and Torque are
        typically body frame, but this is not required.
        
        # Arguments
        * `airstate` - The vehicle airstate at the current timestep
        * `rates` - The body axis rates at the current timestep \[roll,pitch,yaw\] (rad/s)
        * `inputstate` - A reference to the inputstate passed to [AffectedBody::step]
        """
        pass


class AffectedBody:
    """Represent a body subject to aerodynamic effects"""
    def __init__(self, body: AeroBody, effectors: List[AeroEffect]):
        self.body = body # Underlying AeroBody
        self.effectors = effectors # Vec of aerodynamic effects

    def step(self, delta_t: float, inputstate: List[float]) -> None:
        """
        Propagate the system state by delta_t with `inputstate`

        NB: Forces and Torques are calculated at the beginning of the timestep and are not recalculated
        as part of the Runge-Kutta iteration.

        # Arguments
        * `delta_t` - The timestep for this step
        * `inputstate` - The input state to pass to the suplied [AeroEffect]s
        """
        airstate = self.body.get_airstate()
        rates = self.body.rates
        ft_pairs = [ e.get_effect(airstate,rates,inputstate) for e in self.effectors ]
        
        forces,torques = zip(*ft_pairs)
        
        self.body.step(forces,torques,delta_t)
    
    @property
    def acceleration(self) -> List[float]:
        """
        Get body-frame acceleration at the start of the previous timestep

        See [Body::acceleration] for more details
        """
        self.body.acceleration
    
    def set_state(self,new_state: np.ndarray) -> None:
        """
        Set the statevector for the underlying [Body]
    
        The statevector is in the order: \[position,velocity(body),attitude_quaternion(i,j,k,w),axis_rates(body)\]
        
        The statevector must be a column vector
        """
        self.body.set_state(new_state)

    # Implement properties from stateview
    @property
    def position(self) -> List[float]:
        """Return the world-frame position"""
        return self.body.position
    
    @property
    def velocity(self) -> List[float]:
        """Return the body-frame velocity"""
        return self.body.velocity
    
    @property
    def attitude(self) -> List[float]:
        """Return the attitude quaternion"""
        return self.body.attitude
    
    @property
    def rates(self) -> List[float]:
        """Return the body-frame axis rates"""
        return self.body.rates
   
    @property
    def statevector(self) -> List[float]:
        """Return the full statevector"""
        return self.body.statevector

    @property
    def airstate(self) -> AirState:
        """Return the airstate"""
        return self.body.get_airstate()
