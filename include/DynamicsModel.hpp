#ifndef _DYNAMICSMODEL_H
#define _DYNAMICSMODEL_H



// This is meant as a template example of a DynamicsModel
// that is required by the filtering algorithms. You can either
// inherit from this class, or roll your own (probably better speed)
// 
// The required functions are:
//		-- Vector propagate(dt, state, control, process_noise)
//		-- Matrix get_stm(dt, state)
template <typename Vector, typename Matrix>
class DynamicsModel{
public:
	virtual Matrix get_stm(double dt, const Vector & state) const = 0;
	virtual Vector propagate(double dt, const Vector & state, const Vector & control, const Vector & process_noise) const = 0;
};

#endif