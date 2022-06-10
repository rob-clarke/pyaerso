from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import List
import numpy as np


def list_to_column_vec(li: List[float]) -> np.ndarray:
    return np.array([ [e] for e in li ])

class Frame(Enum):
    """Represent reference frame of quantities"""
    WORLD = 0
    BODY = 1


@dataclass
class Force:
    """Represent a force in either body or world reference frame"""
    force: np.ndarray
    frame: Frame
    
    @staticmethod
    def world(x: (float|List[float]),y: (float|None)=None,z: (float|None)=None) -> Force:
        """Create a world-referenced force with components `x`,`y`,`z`"""
        if y is None and z is None:
            return Force.world_vec(x)
        return Force(
            np.array([[x],[y],[z]]),
            Frame.WORLD
        )
        
    @staticmethod
    def world_vec(vec: (List[float]|np.ndarray)) -> Force:
        """Create a world-referenced force with vec"""
        if not isinstance(vec,np.ndarray):
            vec = list_to_column_vec(vec)
        return Force(vec,Frame.WORLD)
        
    @staticmethod
    def body(x: (float|List[float]),y: (float|None)=None,z: (float|None)=None) -> Force:
        """Create a body-referenced force with components `x`,`y`,`z`"""
        if y is None and z is None:
            return Force.body_vec(x)
        return Force(
            np.array([[x],[y],[z]]),
            Frame.BODY
        )
        
    @staticmethod
    def body_vec(vec: (List[float]|np.ndarray)) -> Force:
        """Create a body-referenced force with vec"""
        if not isinstance(vec,np.ndarray):
            vec = list_to_column_vec(vec)
        return Force(vec,Frame.BODY)
        

@dataclass
class Torque:
    """Represent a torque in either body or world reference frame"""
    torque: np.ndarray
    frame: Frame
    
    @staticmethod
    def world(x: (float|List[float]),y: (float|None)=None,z: (float|None)=None) -> Torque:
        """Create a world-referenced torque with components `x`,`y`,`z`"""
        if y is None and z is None:
            return Torque.world_vec(x)
        return Torque(
            np.array([[x],[y],[z]]),
            Frame.WORLD
        )
        
    @staticmethod
    def world_vec(vec: (List[float]|np.ndarray)) -> Torque:
        """Create a world-referenced torque with vec"""
        if not isinstance(vec,np.ndarray):
            vec = list_to_column_vec(vec)
        return Torque(vec,Frame.WORLD)
        
    @staticmethod
    def body(x: (float|List[float]),y: (float|None)=None,z: (float|None)=None) -> Torque:
        """Create a body-referenced torque with components `x`,`y`,`z`"""
        if y is None and z is None:
            return Torque.body_vec(x)
        return Torque(
            np.array([[x],[y],[z]]),
            Frame.BODY
        )
        
    @staticmethod
    def body_vec(vec: (List[float]|np.ndarray)) -> Torque:
        """Create a body-referenced torque with vec"""
        if not isinstance(vec,np.ndarray):
            vec = list_to_column_vec(vec)
        return Torque(vec,Frame.BODY)


class StateView:
    """Class to allow more meaningful access to the state vector"""
    def __init__(self,state: np.ndarray):
        self.state = state
    
    @property
    def position(self) -> np.ndarray:
        """Return the world-frame position"""
        return self.state[0:3,:]
    
    @property
    def velocity(self) -> np.ndarray:
        """Return the body-frame velocity"""
        return self.state[3:6,:]
    
    @property
    def attitude(self) -> np.ndarray:
        """Return the attitude quaternion"""
        return self.state[6:10,:]
    
    @property
    def rates(self) -> np.ndarray:
        """Return the body-frame axis rates"""
        return self.state[10:13,:]

    def __getitem__(self,slice: (int|slice)):
        """Allow access to underlying vector"""
        return self.state.__getitem__(slice)
