//
//  Flow.cpp
//  i-eel
//
//  Created by Jeremie Bec on 25/11/2020.
//

#include "Flow.hpp"

Flow::Flow(double sig) {
	type_ = Shear;
	shear_ = new double[9];
	for(int idim=0; idim<9; ++idim)
	shear_[idim] = 0;
	shear_[1] = sig; // \partial_y u_x = sig
}

double Flow::velocity(double x, double y, double z, int dim) {
	double v = 0.0;
	assert(dim>=0 && dim<3);
	switch (type_) {
		case Null:
			break;
		case Shear:
			v = shear_[3*dim]*x+shear_[3*dim+1]*y+shear_[3*dim+2]*z;
			break;
		default:
			ErrorMsg("Flow type not found");
			break;
	}
	return v;
}

double Flow::gradient(double x, double y, double z, int dim) {
	double g = 0.0;
	assert(dim>=0 && dim<9);
	switch (type_) {
		case Null:
			break;
		case Shear:
			g = shear_[dim];
			break;
		default:
			ErrorMsg("Flow type not found");
			break;
	}
	return g;
}
