#[cfg(not(feature="single-precision"))]
/// Default float representation used in the library
/// 
/// Used when "single-precision" feature is not enabled
pub(crate) type DefaultFloatRepr = f64;

#[cfg(feature="single-precision")]
/// Default float representation used in the library
/// 
/// Used when "single-precision" feature is enabled
pub(crate) type DefaultFloatRepr = f32;
