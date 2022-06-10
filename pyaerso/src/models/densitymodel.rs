use pyo3::prelude::*;

use crate::types::DefaultFloatRepr as FpR;

use aerso::{
    DensityModel,
    types::Vector3
};

#[pyclass(name="DensityModel")]
pub(crate) struct PyDensityModel {
    model: Py<PyAny>
}

impl DensityModel for PyDensityModel {
    fn get_density(&self, position: &Vector3) -> FpR {
        let position_py = [position[0],position[1],position[2]];
        Python::with_gil(|py| {
            let pyresult = self.model.call_method1(py,"get_density",(position_py,));
            if let Err(error) = pyresult {
                panic!("Error calling `get_density`: {}",error);
            }
            if let Ok(result) = pyresult.unwrap().extract(py) {
                result
                }
            else {
                panic!("`get_density` should return a float");
                }
        })
    }
}

impl<'source> FromPyObject<'source> for PyDensityModel {
    fn extract(ob: &'source PyAny) -> Result<PyDensityModel,PyErr> {
        if !ob.hasattr("get_density")? {
            Err(pyo3::exceptions::PyTypeError::new_err("DensityModel needs `get_density` method"))
        }
        else {
            Ok(PyDensityModel { model: ob.into() })
        }
    }
}
