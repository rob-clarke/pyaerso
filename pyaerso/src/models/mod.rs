use crate::types::DefaultFloatRepr as FpR;

mod windmodel;
mod densitymodel;

pub(crate) use windmodel::PyWindModel;
pub(crate) use densitymodel::PyDensityModel;

use aerso::{
    WindModel,DensityModel,
    types::Vector3
};

trait ModelFacade<M,F> {
    fn new(model: M) -> F;
}

pub(crate) struct WindModelFacade {
    pub(crate) model: Box<dyn WindModel>
}

impl WindModel for WindModelFacade {
    fn get_wind(&self, position: &Vector3) -> Vector3 {
        self.model.get_wind(position)
    }
    fn step(&mut self, delta_t: FpR) {
        self.model.step(delta_t)
    }
}

impl ModelFacade<Box<dyn WindModel>,WindModelFacade> for WindModelFacade {
    fn new(model: Box<dyn WindModel>) -> Self {
        Self {
            model
        }
    }
}

pub(crate) struct DensityModelFacade {
    pub(crate) model: Box<dyn DensityModel>
}

impl DensityModel for DensityModelFacade {
    fn get_density(&self, position: &Vector3) -> FpR {
        self.model.get_density(position)
    }
}

impl ModelFacade<Box<dyn DensityModel>,DensityModelFacade> for DensityModelFacade {
    fn new(model: Box<dyn DensityModel>) -> Self {
        Self {
            model
        }
    }
}
