# NPQ_analysis
NPQ relaxation kinetics analysis

Takes Chlorophyll fluorescence time series data and returns non-photochemical quenching relaxation parameters

NPQformat.m formats input data for NPQ main model
NPQmainmodel.m fits a 3-component multi-exponential curve to NPQ. requires either lsqnonlin or lsqcurvefit (matlab optimization package)


u2param.m maps optimized components of model fit to recognizable NPQ relaxation parameters
param2u.m inverse fxn of u2param, maps NPQ relaxation parameters to an optimization friendly vector, u.


