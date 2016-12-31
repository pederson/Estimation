#ifndef _CONTROLMODEL_H
#define _CONTROLMODEL_H

#include "../LinearAlgebra/Matrix.hpp"

class ControlModel{
public:
	virtual Vector get_control(double dt, const Vector & state) const = 0;
	virtual Matrix get_gain(double dt, const Vector & state) const = 0;
};

#endif