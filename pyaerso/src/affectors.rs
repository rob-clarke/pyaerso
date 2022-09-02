use pyo3::prelude::*;

use aerso::{AeroEffect,AirState};
use aerso::types::{Vector3,Force,Torque,StateView};

use crate::types::DefaultFloatRepr as FpR;
use crate::{PyForce,PyTorque};

#[pyclass(name="AeroEffect")]
pub(crate) struct PyAeroEffect {
    pyobj: Py<PyAny>,
}

impl AeroEffect<Vec<FpR>> for PyAeroEffect {
    fn get_effect(&self, airstate: AirState, rates: Vector3, inputstate: &Vec<FpR>) -> (Force,Torque) {
        let airstate_py = [airstate.alpha,airstate.beta,airstate.airspeed,airstate.q];
        let rates_py = [rates.x,rates.y,rates.z];
        
        Python::with_gil(|py| {
            let inputstate_py = inputstate.clone().into_py(py);
            
            let call_result = self.pyobj.call_method1(py,"get_effect",(airstate_py,rates_py,inputstate_py));
            if let Err(error) = call_result {
                panic!("Error calling `get_effect`: {}",error);
            }
            if let Ok(result) = call_result.unwrap().extract(py) {
                let py_forcetorque: (PyRef<PyForce>,PyRef<PyTorque>) = result;
                let force: Force = py_forcetorque.0.clone().into();
                let torque: Torque = py_forcetorque.1.clone().into();
                (force,torque)
            }
            else {
                panic!("`get_effect` should return a (Force,Torque) tuple");
            }
        })
    }
}

impl<'source> FromPyObject<'source> for PyAeroEffect {
    fn extract(ob: &'source PyAny) -> PyResult<PyAeroEffect> {
        if !ob.hasattr("get_effect")? {
            Err(pyo3::exceptions::PyTypeError::new_err("AeroEffect needs `get_effect` method"))
        }
        else {
            Ok(PyAeroEffect { pyobj: ob.into() })
        }
    }
}

use crate::{
    PyAeroBody,
    models::{WindModelFacade,DensityModelFacade}
};

#[pyclass(name="AffectedBody",unsendable)]
pub(crate) struct PyAffectedBody {
    affectedbody: aerso::AffectedBody<Vec<FpR>,FpR,WindModelFacade,DensityModelFacade>,
}

#[pymethods]
impl PyAffectedBody {
    #[new]
    fn new(aerobody: &PyCell<PyAeroBody>, effectors: Vec<PyAeroEffect>) -> PyResult<Self> {
        let aerobody = aerobody.replace(PyAeroBody { aerobody: None });
        Ok(PyAffectedBody {
            affectedbody: aerso::AffectedBody {
                body: aerobody.aerobody.unwrap(),
                effectors: effectors.into_iter().map(|v| Box::new(v) as Box<dyn aerso::AeroEffect>).collect()
            }
        })
    }
    
    pub fn get_derivative(&self, py_state: [FpR;13], input: Vec<FpR>) -> PyResult<[FpR;13]>{
        let state = aerso::types::StateVector::from_vec(Vec::from(py_state));
        Ok(self.affectedbody.get_derivative(&state, &input).into())
    }
    
    pub fn step(&mut self, delta_t: FpR, input: Vec<FpR>) {
        self.affectedbody.step(delta_t, &input)
    }
    
    #[getter]
    fn get_position(&self) -> PyResult<[FpR;3]> {
        Ok(self.affectedbody.position().into())
    }
    
    #[getter]
    fn get_velocity(&self) -> PyResult<[FpR;3]> {
        Ok(self.affectedbody.velocity().into())
    }
    
    #[getter]
    fn get_attitude(&self) -> PyResult<[FpR;4]> {
        let q = self.affectedbody.attitude();
        Ok([q.i, q.j, q.k, q.w])
    }
    
    #[getter]
    fn get_rates(&self) -> PyResult<[FpR;3]> {
        Ok(self.affectedbody.rates().into())
    }
    
    #[getter]
    fn get_statevector(&self) -> PyResult<[FpR;13]> {
        Ok(self.affectedbody.statevector().into())
    }
    
    #[setter]
    fn set_statevector(&mut self, state: [FpR;13]) -> PyResult<()> {
        self.affectedbody.set_state(state.into());
        Ok(())
    }
    
    #[getter]
    fn get_airstate(&self) -> PyResult<[FpR;4]> {
        let airstate = self.affectedbody.get_airstate();
        Ok([airstate.alpha, airstate.beta, airstate.airspeed, airstate.q])
    }
    
}