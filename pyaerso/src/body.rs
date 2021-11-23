use pyo3::prelude::*;

use aerso::types::StateView;

use crate::{PyForce,PyTorque};

#[pyclass(name="Body")]
pub struct PyBody {
    pub(crate) body: aerso::Body
    }

#[pymethods]
impl PyBody {
    
    #[getter]
    fn get_position(&self) -> PyResult<[f64;3]> {
        Ok(self.body.position().into())
    }
    
    #[getter]
    fn get_velocity(&self) -> PyResult<[f64;3]> {
        Ok(self.body.velocity().into())
    }
    
    #[getter]
    fn get_attitude(&self) -> PyResult<[f64;4]> {
        let q = self.body.attitude();
        Ok([q.i, q.j, q.k, q.w])
    }
    
    #[getter]
    fn get_rates(&self) -> PyResult<[f64;3]> {
        Ok(self.body.rates().into())
    }
    
    #[getter]
    fn get_statevector(&self) -> PyResult<[f64;13]> {
        Ok(self.body.statevector().into())
    }
    
    #[new]
    fn new(mass: f64, inertia_py: Vec<Vec<f64>>, position_py: Vec<f64>, velocity_py: Vec<f64>, attitude_py: Vec<f64>, rates_py: Vec<f64>) -> Self {
        let inertia = aerso::types::Matrix3::new(
            inertia_py[0][0],inertia_py[1][0],inertia_py[2][0],
            inertia_py[0][1],inertia_py[1][1],inertia_py[2][1],
            inertia_py[0][2],inertia_py[1][2],inertia_py[2][2]
            );
        
        let statevector = aerso::types::StateVector::from_vec(vec![
            position_py[0], position_py[1], position_py[2],
            velocity_py[0], velocity_py[1], velocity_py[2],
            attitude_py[3], attitude_py[0], attitude_py[1], attitude_py[2],
            rates_py[0],    rates_py[1],    rates_py[2]
            ]);
        
        Self {
            body: aerso::Body::new_from_statevector(mass,inertia,statevector)
        }
    }
    
    fn step(&mut self, forces_py: Vec<PyRef<PyForce>>, torques_py: Vec<PyRef<PyTorque>>, delta_t: f64) {
        let (forces,torques) = crate::force_torque::convert_force_torque(forces_py, torques_py);        
        self.body.step(&forces[..], &torques[..], delta_t);
    }
   
}
