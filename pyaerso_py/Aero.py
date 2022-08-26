from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np

from typing import List, Union

from .types import StateView, Torque, Force, list_to_column_vec
from .Body import Body

class WindModel(ABC):
    """Trait for general wind model"""
    
    @abstractmethod
    def get_wind(self, position: List[float]) -> List[float]:
        """
        Return the current wind at the specified position in world frame coordinates
        Both vectors should be in North-East-Down frame
        """
        pass
    
    def step(self, delta_t: float) -> None:
        """Advance time of the wind model by `delta_t` seconds"""
        pass


class ConstantWind(WindModel):
    """Built-in [WindModel] to represent a constant wind"""
    
    def __init__(self,wind: List[float]):
        """
        Create a new ConstantWind model with `wind`
        
        # Arguments
        
        * `wind` - The wind direction vector (N,E,D)
        """
        self.wind = wind

    def get_wind(self, position: List[float]) -> List[float]:
        return self.wind
    
    def step(self, delta_t: float) -> None:
        pass


class LogWind(WindModel):
    """
    Built-in [WindModel] to represent a log wind profile
    
    https://en.wikipedia.org/wiki/Log_wind_profile
    """
    
    def __init__(self, d: float, z0: float, u_star: float, bearing: float):
        """
        Create a new LogWind with specified parameters
        
        # Arguments
        
        * `d` - Zero plane displacement (m)
        * `z0` - Surface roughness (m)
        * `u_star` - Friction velocity (m·s<sup>-1</sup>)
        * `bearing` - The bearing for the calculated wind vector (deg)
        """
        self.d = d
        self.z0 = z0
        self.u_star = u_star
        self.bearing = bearing

    def get_wind(self, position: List[float]) -> List[float]:
        k = 0.41
        velocity = self.u_star/k * math.ln((position[2] - self.d) / self.z0)
        bearing_rad = math.radians(self.bearing)
        return [
            velocity * math.cos(bearing_rad),
            velocity * math.sin(bearing_rad),
            0.0
        ]
    
    def step(self, delta_t: float) -> None:
        pass


class PowerWind(WindModel):
    """
    Built-in [WindModel] to represent a wind profile power law
    
    https://en.wikipedia.org/wiki/Wind_profile_power_law
    """
    
    def __init__(self, u_r: float, z_r: float, bearing: float, alpha: float):
        """
        Create a new [PowerWind] model with specified parameters

        # Arguments

        * `u_r` - Reference wind speed (m·s<sup>-1</sup>) 
        * `z_r` - Reference wind height (m)
        * `bearing` - The bearing for the calculated wind vector (deg)
        * `alpha` - Power law exponent (Typically 0.143)
        """
        self.u_r = u_r
        self.z_r = z_r
        self.bearing = bearing
        self.alpha = alpha

    def get_wind(self, position: List[float]) -> List[float]:
        velocity = self.u_r * math.pow(position[2] / self.z, self.alpha)
        bearing_rad = math.radians(self.bearing)
        return [
            velocity * math.cos(bearing_rad),
            velocity * math.sin(bearing_rad),
            0.0
        ]
    
    def step(self, delta_t: float) -> None:
        pass

class DensityModel(ABC):
    """Trait for general density model"""
    
    def get_density(self, position: List[float]) -> float:
        """Return the current density at the specified position (kg.m^-3)"""
        pass


class ConstantDensity(DensityModel):
    """
    Built-in [DensityModel] for constant ISA standard sea level density

    This model does not vary density with altitude.
    """

    def get_density(self, position: List[float]) -> float:
        return 1.225


class StandardDensity(DensityModel):
    """
    Built-in [DensityModel] for ISA standard density model

    This model is valid up to 11km altitude

    http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/atmos/atmos.html
    """
    T_LR = 0.0065; # K/m
    R    = 287.0;  # m^2/s^2/K
    ISA_SL_T  = 288.15; # K
    ISA_11K_T = 216.65; # K
    ISA_SL_H  =      0.0; # m
    ISA_11K_H = 11_000.0; # m
    ISA_SL_P = 101_325.0; # Pa
    
    @staticmethod
    def _interp(x: float, x_0: float, x_1: float, y_0: float, y_1: float) -> float:
        return y_0 + (x-x_0)/(x_1-x_0)*(y_1-y_0)
    
    @staticmethod
    def _get_standard_temperature(altitude: float) -> float:
        return StandardDensity._interp(
            altitude,
            StandardDensity.ISA_SL_H, StandardDensity.ISA_11K_H,
            StandardDensity.ISA_SL_T, StandardDensity.ISA_11K_T
            )
    
    @staticmethod
    def _get_standard_pressure(altitude: float) -> float:
        t_0 = StandardDensity.ISA_SL_T
        p_0 = StandardDensity.ISA_SL_P
        t_a = StandardDensity._get_standard_temperature(altitude)
        lr = StandardDensity.T_LR
        r_dry = StandardDensity.R
        g = 9.81
        
        return p_0 * (t_a/t_0) ** ( g/(lr*r_dry) )
    
    def get_density(self, position: List[float]) -> float:
        altitude = -position[2]
        r_dry = StandardDensity.R
        return StandardDensity._get_standard_pressure(altitude) / (r_dry * StandardDensity._get_standard_temperature(altitude))


@dataclass
class AirState:
    """Represent generic air state"""
    alpha: float # Angle of attack (radians)
    beta: float # Angle of sideslip (radians)
    airspeed: float # Airspeed (m·s<sup>-1</sup>) 
    q: float # Dynamic pressure (Pa) (kg·m<sup>-1</sup>·s<sup>2</sup>)
    
    def __getitem__(self,slice: Union[int,slice]) -> float:
        return ([self.alpha,self.beta,self.airspeed,self.q])[slice]


class AeroBody:
    """Represent a body in an atmosphere"""
    def __init__(self, body: Body, wind_model: WindModel=None, density_model: DensityModel=None):
        if wind_model is None:
            wind_model = ConstantWind( np.zeros((3,1)) )
        elif isinstance(wind_model,tuple):
            name = wind_model[0]
            args = wind_model[1]
            if name == "ConstantWind":
                wind_model = ConstantWind([*args])
            if name == "LogWind":
                wind_model = LogWind(*args)
            if name == "PowerWind":
                wind_model = PowerWind(*args)
        
        if density_model is None:
            density_model = ConstantDensity()
        elif isinstance(density_model,tuple):
            name = density_model[0]
            args = density_model[1]
            if name == "StandardDensity":
                density_model = StandardDensity(*args)
            if name == "ConstantDensity":
                density_model = ConstantDensity(*args)
        
        self.body = body # The underlying rigid body
        self.wind_model = wind_model # Optional wind model
        self.density_model = density_model # Optional density model

    @staticmethod
    def new(body: Body) -> AeroBody:
        """
        Create an [AeroBody] with no wind and constant ISA standard sea-level density
        
        # Arguments
        
        * `body` - The kinematics body to use
        """
        return AeroBody.with_wind_model(
            body,
            ConstantWind( np.zeros((3,1)) )
            )

    @staticmethod
    def with_wind_model(body: Body, wind_model: WindModel) -> AeroBody:
        """
        Create an AeroBody with a [WindModel] and constant ISA standard sea-level density

        # Arguments

        * `body` - The kinematics body to use
        * `wind_model` - The [WindModel] to use
        """
        density_model = StandardDensity()
        return AeroBody.with_density_model(body,wind_model,density_model)

    @staticmethod
    def with_density_model(body: Body, wind_model: WindModel=None, density_model: DensityModel=None):
        """
        Create an AeroBody with a [WindModel] and a [DensityModel]

        # Arguments

        * `body` - The kinematics body to use
        * `wind_model` - The [WindModel] to use
        * `density_model` - The [DensityModel] to use
        """
        return AeroBody(body,wind_model,density_model)

    def get_airstate(self) -> AirState:   
        """
        Return an [AirState] representing the current aerodynamic state of the body

        The [AirState] includes the angles of attack (`alpha`) and sideslip (`beta`), the `airspeed` and the dynamic pressure, (`q`).

        It is calculated using the supplied wind and density models.
        """
        
        current_world_wind = self.wind_model.get_wind(self.body.position)
        current_body_wind = list_to_column_vec(self.body.velocity) - (Body.get_dcm(StateView(self.body.statevector)) @ current_world_wind)
        
        [u,v,w] = current_body_wind[:,0]
        
        u_sqd = u**2
        v_sqd = v**2
        w_sqd = w**2
        
        airspeed = np.sqrt( u_sqd + v_sqd + w_sqd )
        
        alpha = np.arctan2(w,u)
        
        beta = np.arcsin(v/airspeed) if airspeed > 0.001 else 0.0
        
        q = 0.5 * self.density_model.get_density(self.body.position) * airspeed**2
        
        return AirState(alpha, beta, airspeed, q )
    
    def step(self, forces: List[Force], torques: List[Torque], delta_t: float) -> None:
        """
        Propagate the body state and wind_model by `delta_t` under the supplied `forces` and `torques`
        
        See the documentation for [Body::step] for further details
        """
        self.body.step(forces, torques, delta_t)
        self.wind_model.step(delta_t)
    
    
    @property
    def acceleration(self) -> List[float]:
        """
        Get body-frame acceleration at the start of the previous timestep

        See [Body::acceleration] for more details
        """
        self.body.acceleration
    
    def set_state(self,new_state) -> None:
        """
        Set the statevector for the underlying [Body]
    
        The statevector is in the order: \[position,velocity(body),attitude_quaternion(i,j,k,w),axis_rates(body)\]
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
        return self.body.statevector[:,0]
