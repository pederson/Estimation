#ifndef _PROCESSNOISEMODEL_H
#define _PROCESSNOISEMODEL_H


// This is meant as a template example of a ProcessNoiseModel
// that is required by the filtering algorithms. You can either
// inherit from this class, or roll your own (probably better speed)
// 
// The required functions are:
//		-- Vector get_mean(dt, state)
//		-- Matrix get_pnmt(dt, state)
//		-- Matrix get_covariance(dt, state)
//
// Note that not all functions are required for every filter
template <typename Vector, typename Matrix>
class ProcessNoiseModel{
public:
	virtual Vector get_mean(double dt, const Vector & state) const = 0;
	virtual Matrix get_covariance(double dt, const Vector & state) const = 0;
	virtual Matrix get_pntm(double dt, const Vector & state) const = 0;
};

#endif