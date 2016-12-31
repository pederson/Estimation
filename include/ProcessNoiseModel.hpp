#ifndef _PROCESSNOISEMODEL_H
#define _PROCESSNOISEMODEL_H

#include "../LinearAlgebra/Matrix.hpp"

class ProcessNoiseModel{
public:
	virtual Vector get_mean(double dt, const Vector & state) const = 0;
	virtual Matrix get_covariance(double dt, const Vector & state) const = 0;
	virtual Matrix get_pntm(double dt, const Vector & state) const = 0;
};

#endif