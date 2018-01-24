#ifndef _CONTROLMODEL_H
#define _CONTROLMODEL_H


// This is meant as a template example of a ControlModel
// that is required by the filtering algorithms. You can either
// inherit from this class, or roll your own (probably better speed)
// 
// The required functions are:
//		-- Vector get_control(dt, state)
//		-- Matrix get_gain(dt, state)
template <typename Vector, typename Matrix>
class ControlModel{
public:
	virtual Vector get_control(double dt, const Vector & state) const = 0;
	virtual Matrix get_gain(double dt, const Vector & state) const = 0;
};

#endif