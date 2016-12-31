#ifndef _MEASUREMENTMODEL_H
#define _MEASUREMENTMODEL_H

#include "../LinearAlgebra/Matrix.hpp"

class MeasurementModel{
public:
	virtual Matrix get_jacobian(double dt, const Vector & state) const = 0;
	virtual Matrix get_covariance(double dt, const Vector & state) const = 0;
};

#endif