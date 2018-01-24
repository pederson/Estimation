#ifndef _MEASUREMENTMODEL_H
#define _MEASUREMENTMODEL_H


// This is meant as a template example of a MeasurementModel
// that is required by the filtering algorithms. You can either
// inherit from this class, or roll your own (probably better speed)
// 
// The required functions are:
//		-- Vector get_value(dt, state)
//		-- Matrix get_jacobian(dt, state)
//		-- Matrix get_covariance(dt, state)
//		-- Vector get_residual(state, measurement)
//
// Note that not all functions are required for every filter
template <typename Vector, typename Matrix>
class MeasurementModel{
public:
	virtual Vector get_value(double dt, const Vector & state) const = 0;
	virtual Matrix get_jacobian(double dt, const Vector & state) const = 0;
	virtual Matrix get_covariance(double dt, const Vector & state) const = 0;
	virtual Vector get_residual(const Vector & state, const Vector & measurement) const = 0;
};

#endif