use pyo3::prelude::*;

use aerso::types::StateView;

use crate::{PyForce,PyTorque};
use crate::types::DefaultFloatRepr as FpR;

#[pyclass(name="Body")]
pub struct PyBody {
    pub(crate) body: aerso::Body
    }

#[pymethods]
impl PyBody {
    
    #[getter]
    fn get_position(&self) -> PyResult<[FpR;3]> {
        Ok(self.body.position().into())
    }
    
    #[getter]
    fn get_velocity(&self) -> PyResult<[FpR;3]> {
        Ok(self.body.velocity().into())
    }
    
    #[getter]
    fn get_attitude(&self) -> PyResult<[FpR;4]> {
        let q = self.body.attitude();
        Ok([q.i, q.j, q.k, q.w])
    }
    
    #[getter]
    fn get_rates(&self) -> PyResult<[FpR;3]> {
        Ok(self.body.rates().into())
    }
    
    #[getter]
    fn get_statevector(&self) -> PyResult<[FpR;13]> {
        Ok(self.body.statevector().into())
    }
    
    #[setter]
    fn set_statevector(&mut self, state: [FpR;13]) -> PyResult<()> {
        self.body.set_state(state.into());
        Ok(())
    }
    
    #[new]
    fn new(mass: FpR, inertia_py: Vec<Vec<FpR>>, position_py: Vec<FpR>, velocity_py: Vec<FpR>, attitude_py: Vec<FpR>, rates_py: Vec<FpR>) -> Self {
        let inertia = aerso::types::Matrix3::new(
            inertia_py[0][0],inertia_py[1][0],inertia_py[2][0],
            inertia_py[0][1],inertia_py[1][1],inertia_py[2][1],
            inertia_py[0][2],inertia_py[1][2],inertia_py[2][2]
            );
        
        let statevector = aerso::types::StateVector::from_vec(vec![
            position_py[0], position_py[1], position_py[2],
            velocity_py[0], velocity_py[1], velocity_py[2],
            attitude_py[0], attitude_py[1], attitude_py[2], attitude_py[3],
            rates_py[0],    rates_py[1],    rates_py[2]
            ]);
        
        Self {
            body: aerso::Body::new_from_statevector(mass,inertia,statevector)
        }
    }

    fn step(&mut self, forces_py: Vec<PyRef<PyForce>>, torques_py: Vec<PyRef<PyTorque>>, delta_t: FpR) {
        let (forces,torques) = crate::force_torque::convert_force_torque(forces_py, torques_py);        
        self.body.step(&forces[..], &torques[..], delta_t);
    }

    fn __str__(&self) -> String {
        String::from(
            format!("Body()")
        )
    }

    fn __repr__(&self) -> String {
        String::from(
            format!("Body()")
        )
    }
}
