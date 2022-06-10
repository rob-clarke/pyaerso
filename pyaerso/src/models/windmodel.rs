use pyo3::prelude::*;

use crate::types::DefaultFloatRepr as FpR;

use aerso::{
    WindModel,
    types::Vector3
};

#[pyclass(name="WindModel")]
pub(crate) struct PyWindModel {
    model: Py<PyAny>,
}

impl WindModel for PyWindModel {
    fn get_wind(&self, position: &Vector3) -> Vector3 {
        let position_py = [position[0],position[1],position[2]];
        let pyvec: [FpR;3] = Python::with_gil(|py| {
            let pyresult = self.model.call_method1(py,"get_wind",(position_py,));
            if let Err(error) = pyresult {
                panic!("Error calling `get_wind`: {}",error);
            }
            if let Ok(result) = pyresult.unwrap().extract(py) {
                result
                }
            else {
                panic!("`get_wind` should return a List of 3 floats for N,E,D components of the wind vector");
                }
        });

        Vector3::new(pyvec[0],pyvec[1],pyvec[2])
    }
    
    fn step(&mut self, delta_t: FpR) {
        Python::with_gil(|py| {
            let pyresult = self.model.call_method1(py,"step",(delta_t,));
            if let Err(error) = pyresult {
                panic!("Error calling `step`: {}",error);
            }
        });
    }
}

impl<'source> FromPyObject<'source> for PyWindModel {
    fn extract(ob: &'source PyAny) -> Result<PyWindModel,PyErr> {
        if !ob.hasattr("get_wind")? || !ob.hasattr("step")? {
            Err(pyo3::exceptions::PyTypeError::new_err("WindModel needs methods `get_wind` and `step`"))
        }
        else {
            Ok(PyWindModel { model: ob.into() })
        }
    }
}
