use std::ops::Deref;

use pyo3::prelude::*;
use pyo3::types::PyLong;
use pyo3::exceptions;

use aerso::types::{Force,Torque};
use crate::types::DefaultFloatRepr as FpR;

#[derive(Copy,Clone,PartialEq)]
pub enum PyFrame {
    World = aerso::types::Frame::World as isize,
    Body = aerso::types::Frame::Body as isize,
}

impl<'source> FromPyObject<'source> for PyFrame {
    fn extract(ob: &'source PyAny) -> Result<PyFrame,PyErr> {
        match ob.cast_as::<PyLong>()?.extract::<isize>()? {
            x if x == PyFrame::World as isize => Ok(PyFrame::World),
            x if x == PyFrame::Body as isize => Ok(PyFrame::Body),
            _ => Err(exceptions::PyTypeError::new_err("Unable to extract valid enum value")),
        }
    }
}

#[pyclass(name="Force")]
#[derive(Clone)]
pub struct PyForce {
    pub force: Vec<FpR>,
    pub frame: PyFrame,
}

impl PyForce {
    fn new(force_py: Vec<FpR>, frame_py: PyFrame) -> Self {        
        PyForce {
            force: vec![force_py[0],force_py[1],force_py[2]],
            frame: frame_py
        }
    }
}

#[pymethods]
impl PyForce {
    
    #[staticmethod]
    fn world(force_py: Vec<FpR>) -> Self {
        PyForce::new(force_py,PyFrame::World)
    }
    
    #[staticmethod]
    fn body(force_py: Vec<FpR>) -> Self {
        PyForce::new(force_py,PyFrame::Body)
    }
    
    #[getter]
    fn force(&self) -> Vec<FpR> {
        vec![self.force[0],self.force[1],self.force[2]]
    }
    
    #[getter]
    fn frame(&self) -> isize {
        self.frame as isize
    }
    
}

impl Into<Force> for PyForce {
    fn into(self) -> Force {
        match self {
            PyForce { force, frame } if frame == PyFrame::Body => Force::body(force[0],force[1],force[2]),
            PyForce { force, frame } if frame == PyFrame::World => Force::world(force[0],force[1],force[2]),
            _ => { unreachable!(); }
        }
    }
}

#[pyclass(name="Torque")]
#[derive(Clone)]
pub struct PyTorque {
    pub torque: Vec<FpR>,
    pub frame: PyFrame,
}

impl PyTorque {
    fn new(torque_py: Vec<FpR>, frame_py: PyFrame) -> Self {        
        PyTorque {
            torque: vec![torque_py[0],torque_py[1],torque_py[2]],
            frame: frame_py,
        }
    }
}

#[pymethods]
impl PyTorque {
    
    #[staticmethod]
    fn world(force_py: Vec<FpR>) -> Self {
        PyTorque::new(force_py,PyFrame::World)
    }
    
    #[staticmethod]
    fn body(force_py: Vec<FpR>) -> Self {
        PyTorque::new(force_py,PyFrame::Body)
    }
    
    #[getter]
    fn torque(&self) -> Vec<FpR> {
        vec![self.torque[0],self.torque[1],self.torque[2]]
    }
    
    #[getter]
    fn frame(&self) -> isize {
        self.frame as isize
    }
}

impl Into<Torque> for PyTorque {
    fn into(self) -> Torque {
        match self {
            PyTorque { torque, frame } if frame == PyFrame::Body => Torque::body(torque[0],torque[1],torque[2]),
            PyTorque { torque, frame } if frame == PyFrame::World => Torque::world(torque[0],torque[1],torque[2]),
            _ => { unreachable!(); }
        }
    }
}

pub fn convert_force_torque(forces_py: Vec<PyRef<PyForce>>, torques_py: Vec<PyRef<PyTorque>>) -> (Vec<Force>,Vec<Torque>) {
    let forces : Vec<Force> = forces_py.iter().map(|v|
        match v.deref() {
            PyForce { force, frame } if *frame == PyFrame::Body => Force::body(force[0],force[1],force[2]),
            PyForce { force, frame } if *frame == PyFrame::World => Force::world(force[0],force[1],force[2]),
            _ => { unreachable!(); }
        }).collect();

    let torques: Vec<Torque> = torques_py.iter().map(|v|
        match v.deref() {
            PyTorque { torque, frame } if *frame == PyFrame::Body => Torque::body(torque[0],torque[1],torque[2]),
            PyTorque { torque, frame } if *frame == PyFrame::World => Torque::world(torque[0],torque[1],torque[2]),
            _ => { unreachable!(); }
        }).collect();
    
    (forces,torques)
}
