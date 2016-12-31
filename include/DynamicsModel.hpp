#ifndef _DYNAMICSMODEL_H
#define _DYNAMICSMODEL_H

#include "../LinearAlgebra/Matrix.hpp"

class DynamicsModel{
public:
	virtual Matrix get_stm(double dt, const Vector & state) const = 0;
	virtual Vector propagate(double dt, const Vector & state, const Vector & control, const Vector & process_noise) const = 0;
};

#endif