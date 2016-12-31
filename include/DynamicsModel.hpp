#ifndef _DYNAMICSMODEL_H
#define _DYNAMICSMODEL_H

#include "../LinearAlgebra/Matrix.hpp"

class DynamicsModel{
public:
	virtual Matrix get_stm(double dt, const Vector & state) const = 0;
};

#endif