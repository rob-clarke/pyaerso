use pyo3::prelude::*;

use crate::PyBody;

use aerso::types::{
    Vector3,
    StateView
};

use crate::{PyForce,PyTorque};

use crate::models::{
    PyWindModel,WindModelFacade,
    PyDensityModel,DensityModelFacade
    };

#[pyclass(name="AeroBody",unsendable)]
pub struct PyAeroBody {
    pub(crate) aerobody: aerso::AeroBody<f64,WindModelFacade,DensityModelFacade>,
}


#[pymethods]
impl PyAeroBody {
    
    /// Create a new (Python) AeroBody
    /// body: Ref to a Python pyaerso.Body
    /// wind_model: A (string,[*args]) tuple naming a default wind model from aerso, or a Python object with `get_wind` and `step` methods
    /// density_model: A (string,[*args]) tuple naming a default density model from aerso, or a Python object with a `get_density` method
    #[new]
    fn new(body: &PyCell<PyBody>, py_wind_model: Option<Py<PyAny>>, py_density_model: Option<Py<PyAny>>) -> PyResult<PyAeroBody> {
        
        // fn get_named_model<P: pyo3::PyClass + Clone, F: ModelFacade<Box<P>,F>>(pyobject: Py<PyAny>) -> PyResult<F> {
        //     if let Ok(modelname) = Python::with_gil(|py| pyobject.extract::<String>(py)) {
        //         return match modelname.as_str() {
        //             _ => { return Err(pyo3::exceptions::PyValueError::new_err("Unknown model name")); },
        //         }
        //     }
        //     let result = Python::with_gil(|py| pyobject.extract::<P>(py))?;
        //     Ok(F::new(Box::new(result)))
        // }
        
        fn get_windmodel(pyobject: Py<PyAny>) -> PyResult<WindModelFacade> {
            if let Ok((modelname,margs)) = Python::with_gil(|py| pyobject.extract::<(String,Vec<f64>)>(py)) {
                return match modelname.as_str() {
                    "ConstantWind" => {
                        if margs.len() < 3 { Err(pyo3::exceptions::PyValueError::new_err("ConstantWind takes 3 arguments")) }
                        else { Ok(WindModelFacade { model: Box::new(aerso::wind_models::ConstantWind::new(Vector3::new(margs[0],margs[1],margs[2]))) }) }
                    },
                    "PowerWind" => {
                        if margs.len() < 3 { Err(pyo3::exceptions::PyValueError::new_err("PowerWind takes 3 arguments")) }
                        else { Ok(WindModelFacade { model: Box::new(aerso::wind_models::PowerWind::new(margs[0],margs[1],margs[2])) }) }
                    },
                    "LogWind" => {
                        if margs.len() < 4 { Err(pyo3::exceptions::PyValueError::new_err("LogWind takes 4 arguments")) }
                        else { Ok(WindModelFacade { model: Box::new(aerso::wind_models::LogWind::new(margs[0],margs[1],margs[2],margs[3])) }) }
                    },
                    _ => Err(pyo3::exceptions::PyValueError::new_err("Unknown wind model name")),
                }
            }
            let result = Python::with_gil(|py| pyobject.extract::<PyWindModel>(py))?;
            Ok(WindModelFacade { model: Box::new(result) })
        }
        
        fn get_densitymodel(pyobject: Py<PyAny>) -> PyResult<DensityModelFacade> {
            if let Ok((modelname,_)) = Python::with_gil(|py| pyobject.extract::<(String,Vec<f64>)>(py)) {
                return match modelname.as_str() {
                    "StandardDensity" => Ok(DensityModelFacade { model: Box::new(aerso::density_models::StandardDensity {}) }),
                    _ => Err(pyo3::exceptions::PyValueError::new_err("Unknown density model name")),
                }
            }
            let result = Python::with_gil(|py| pyobject.extract::<PyDensityModel>(py))?;
            Ok(DensityModelFacade { model: Box::new(result) })
        }
        
        let wind_model: WindModelFacade = match py_wind_model {
            Some(pyobject) => get_windmodel(pyobject)?,
            None => WindModelFacade { model: Box::new(aerso::wind_models::ConstantWind::new(Vector3::new(0.0,0.0,0.0)) ) }
        };
        
        let density_model: DensityModelFacade = match py_density_model {
            Some(pyobject) => get_densitymodel(pyobject)?,
            None => DensityModelFacade { model: Box::new(aerso::density_models::StandardDensity {}) },
        };
        
        Ok(PyAeroBody {
            aerobody: aerso::AeroBody::with_density_model(
                (&mut *body.borrow_mut()).body,
                wind_model,
                density_model),
        })
    }
    
    #[getter]
    fn get_position(&self) -> PyResult<[f64;3]> {
        Ok(self.aerobody.position().into())
    }
    
    #[getter]
    fn get_velocity(&self) -> PyResult<[f64;3]> {
        Ok(self.aerobody.velocity().into())
    }
    
    #[getter]
    fn get_attitude(&self) -> PyResult<[f64;4]> {
        let q = self.aerobody.attitude();
        Ok([q.i, q.j, q.k, q.w])
    }
    
    #[getter]
    fn get_rates(&self) -> PyResult<[f64;3]> {
        Ok(self.aerobody.rates().into())
    }
    
    #[getter]
    fn get_statevector(&self) -> PyResult<[f64;13]> {
        Ok(self.aerobody.statevector().into())
    }
    
    #[getter]
    fn get_airstate(&self) -> PyResult<[f64;4]> {
        let airstate = self.aerobody.get_airstate();
        Ok([airstate.alpha, airstate.beta, airstate.airspeed, airstate.q])
    }
    
    fn step(&mut self, forces_py: Vec<PyRef<PyForce>>, torques_py: Vec<PyRef<PyTorque>>, delta_t: f64) {
        let (forces,torques) = crate::force_torque::convert_force_torque(forces_py, torques_py);        
        self.aerobody.step(&forces[..], &torques[..], delta_t);
    }
    
}
