use pyo3::prelude::*;

mod types;

mod force_torque;
use force_torque::{PyForce,PyTorque};

mod body;
use body::PyBody;

mod models;

mod aero;
use aero::PyAeroBody;

mod affectors;
use affectors::{PyAeroEffect,PyAffectedBody};


#[pyclass(name="Frame")]
struct PyFrameEnum;

#[pymethods]
impl PyFrameEnum {
    #[classattr]
    const WORLD: isize = aerso::types::Frame::World as isize;
    
    #[classattr]
    const BODY: isize = aerso::types::Frame::Body as isize;
}

#[pymodule]
fn pyaerso(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyBody>()?;
    m.add_class::<PyForce>()?;
    m.add_class::<PyTorque>()?;
    m.add_class::<PyFrameEnum>()?;
    m.add_class::<PyAeroBody>()?;
    m.add_class::<PyAeroEffect>()?;
    m.add_class::<PyAffectedBody>()?;
    Ok(())
}
