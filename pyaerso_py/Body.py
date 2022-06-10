from __future__ import annotations

from typing import List
import numpy as np

from .types import Force, Frame, StateView, Torque, list_to_column_vec

# Integrating Rotations using Non-Unit Quaternions
# https://par.nsf.gov/servlets/purl/10097724
CONSTRAIN_QNORM_DRIFT = False

import os
USE_NUMPY_CROSSPRODUCT = "NUMPY_XPROD" in os.environ
if USE_NUMPY_CROSSPRODUCT:
    cross_prod = np.cross
else:
    def cross_prod(a,b,axis=None):
        return np.array([
            [a[1,0]*b[2,0] - a[2,0]*b[1,0]],
            [a[2,0]*b[0,0] - a[0,0]*b[2,0]],
            [a[0,0]*b[1,0] - a[1,0]*b[0,0]]
        ])

class Body:
    """Represent a 6DoF body affected by gravity"""
    
    def __init__(self, mass: float, inertia: List[List[float]],
        position: List[float], velocity: List[float], attitude: List[float], rates: List[float]
    ):
        # Mass of body (kg)
        self.mass = mass
        # Inertia matrix
        self.inertia = inertia
        self.inertia_inverse = np.linalg.inv(inertia)
        # 13-dimensional state vector
        # Statevector is formed of \[position,velocity(body),attitude_quaternion(i,j,k,w),axis_rates(body)\]
        self.statevector = list_to_column_vec([*position,*velocity,*attitude,*rates])
        # Body frame acceleration of vehicle during last step
        self.last_acceleration: np.zeros((1,3))

    @staticmethod
    def new(mass: float, inertia: List[List[float]],
        position: List[float], velocity: List[float], attitude: List[float], rates: List[float]
    ) -> Body:
        """Create a new instance of Body with `mass` and `inertia` in specified state"""
        return Body(mass,inertia,position,velocity,attitude,rates)
    
    @staticmethod
    def new_at_origin(mass: float, inertia: List[List[float]]) -> Body:
        """Create a new instance of Body with `mass` and `inertia` at the origin"""
        return Body.new(mass,inertia,[0,0,0],[0,0,0],[0,0,0,1],[0,0,0])
    
    @staticmethod
    def new_from_statevector(mass: float, inertia: List[List[float]], statevector: List[float]) -> Body:
        """
        Create a new instance of Body with `mass` and `inertia` in specified state
        statevector is made of \[position,velocity(body),attitude_quaternion(i,j,k,w),axis_rates(body)\]
        """
        position = statevector[0:3]
        velocity = statevector[3:6]
        attitude = statevector[6:10]
        rates =    statevector[10:13]
        return Body(mass,inertia,position,velocity,attitude,rates)
    
    @staticmethod
    def get_dcm(state: StateView) -> np.ndarray:
        """
        Construct the Direction Cosine Matrix (DCM) from the state attitude
        
        Transforms quantites from the world frame to the body frame
        
        Note that this is not a struct method. Usage:
        ```
        body = Body.new_at_origin(1.0,np.eye(3))
        dcm = Body.get_dcm(body.statevector)
        ```
        
        # Arguments
        
        * `state` - The statevector to calculate the DCM for
        """
        # Don't use attitude here to avoid unessecary square root call
        [q1,q2,q3,q0] = state[6:10,0]
        
        q02 = q0**2 # Real part (w)
        q12 = q1**2 # i
        q22 = q2**2 # j
        q32 = q3**2 # k
        
        # NB: This appears as the transpose of Eq. (13) from the referenced paper
        # This matches the convention of the Stengel notes
        return np.array([
            [q02+q12-q22-q32,   2.0*(q1*q2+q0*q3), 2.0*(q1*q3-q0*q2)],
            [2.0*(q1*q2-q0*q3), q02-q12+q22-q32,   2.0*(q2*q3+q0*q1)],
            [2.0*(q1*q3+q0*q2), 2.0*(q2*q3-q0*q1), q02-q12-q22+q32]
        ]) * 1.0/(q02+q12+q22+q32)
    
    @staticmethod
    def get_dcm_body(state: StateView) -> np.ndarray:
        """
        Construct the inverse DCM

        Transforms quantities from the body frame to the world frame

        Note that this is not a struct method. Usage:
        ```
        body = Body.new_at_origin(1.0,np.eye(3));
        dcm = Body.get_dcm_body(body.statevector);
        ```

        # Arguments

        * `state` - The statevector to calculate the DCM for
        """
        Body.get_dcm(state).transpose()
    
    def get_derivative(self, state: StateView, forces: List[Force], torques: List[Torque]) -> np.ndarray:
        """
        Calculate the state derivative
        
        NB: Gravity is included by default
        
        # Arguments
        * `state` - 13-dimensional state vector to get derivative about
        * `forces` - Vector of applied forces, both world and body frame
        * `torques` - Vector of applied torques, both world and body frame
        """
        gravity_accel = np.array([[0],[0],[9.81]])
        
        world_forces = gravity_accel * self.mass
        body_forces = np.zeros((3,1))
        for force in forces:
            if force.frame == Frame.WORLD:
                world_forces += force.force
            if force.frame == Frame.BODY:
                body_forces += force.force
        
        world_torques = np.zeros((3,1))
        body_torques = np.zeros((3,1))
        for torque in torques:
            if torque.frame == Frame.WORLD:
                world_torques += torque.torque
            if torque.frame == Frame.BODY:
                body_torques += torque.torque
        
        dcm = Body.get_dcm(state)
        dcm_body = dcm.transpose()

        position_dot = dcm_body @ state.velocity
        velocity_dot = cross_prod(state.velocity,state.rates,axis=0) + ( (dcm @ world_forces) + body_forces ) * 1.0/self.mass
        
        [o_x,o_y,o_z] = state.rates[:,0]
        # NB: This matrix appears different from Eq. (21) in the reference due to quaternion ordering
        # The paper uses [w,i,j,k] but nalgebra uses [i,j,k,w]
        # Hence row/column 1 is shifted to row/column 4
        qdot_matrix = np.array([
            [  0.0,  o_z, -o_y,  o_x ],
            [ -o_z,  0.0,  o_x,  o_y ],
            [  o_y, -o_x,  0.0,  o_z ],
            [ -o_x, -o_y, -o_z,  0.0 ]
        ])
        
        # NB: Quaternion does not remain normalised throughout integration
        q = state[6:10,:] # Don't use attitude here to avoid uneccesary square root
        
        if not CONSTRAIN_QNORM_DRIFT:
            attitude_dot = 0.5 * qdot_matrix @ q
        else:
            # Use Eq. (23) from the paper to set the free parameter to constrain the drift of the quaternion norm
            qnorm_sqd = q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2
            k = 0.001
            c = k * (1.0 - qnorm_sqd)
            attitude_dot = 0.5 * qdot_matrix @ q + q * c
        
        rates_dot = self.inertia_inverse @ (dcm @ world_torques + body_torques - cross_prod((self.inertia @ state.rates),state.rates,axis=0) )
        
        # print(position_dot)
        # print(velocity_dot)
        # print(attitude_dot)
        # print(rates_dot)
        
        return np.array([
            *position_dot,
            *velocity_dot,
            *attitude_dot,
            *rates_dot
            ])
        
     
    def step(self, forces: List[Force], torques: List[Torque], delta_t: float) -> None:
        """
        Propagate the state vector by delta_t under the supplied forces and torques

        Uses 4th-order Runge-Kutta integration

        NB: Gravity is included by default

        # Arguments

        * `forces` - Vector of applied forces, both world and body frame
        * `torques` - Vector of applied torques, both world and body frame
        * `delta_t` - Timestep (s)
        """
        k1 = self.get_derivative(StateView(self.statevector),                    forces, torques)
        k2 = self.get_derivative(StateView(self.statevector + k1 * delta_t/2.0), forces, torques)
        k3 = self.get_derivative(StateView(self.statevector + k2 * delta_t/2.0), forces, torques)
        k4 = self.get_derivative(StateView(self.statevector + k3 * delta_t),     forces, torques)
        
        # NB: k1 is a derivative so velocity -> velocity_dot -> acceleration
        self.last_acceleration = StateView(k1).velocity
        self.statevector += (k1 + k2*2.0 + k3*2.0 + k4) * delta_t/6.0
    
    @property
    def acceleration(self) -> List[float]:
        """
        Get body-frame acceleration at the start of the previous timestep

        The resultant acceleration is in body frame and is the coordinate acceleration.
        To turn this into proper acceleration as seen by an accelerometer, the acceleration of the
        frame needs to be added, in this case standard gravity:
        ```
        acc_frame = Body.get_dcm(body.statevector) * np.array([[0],[0],[-9.81]])
        acc_proper = body.acceleration + acc_frame
        ```
        """
        self.last_acceleration[:,0]
    
    def set_state(self,new_state: np.ndarray) -> None:
        """
        Set the body statevector
        
        Will also reset the body acceleration to zero
        
        The statevector is in the order: \[position,velocity(body),attitude_quaternion(i,j,k,w),axis_rates(body)\]
        
        The statevector must be a column vector
        """
        self.statevector = new_state
        self.last_acceleration = np.zeros((3,1))

    # Implement properties from stateview
    @property
    def position(self) -> List[float]:
        """Return the world-frame position"""
        return self.statevector[0:3,0]
    
    @property
    def velocity(self) -> List[float]:
        """Return the body-frame velocity"""
        return self.statevector[3:6,0]
    
    @property
    def attitude(self) -> List[float]:
        """Return the attitude quaternion"""
        return self.statevector[6:10,0]
    
    @property
    def rates(self) -> List[float]:
        """Return the body-frame axis rates"""
        return self.statevector[10:13,0]
