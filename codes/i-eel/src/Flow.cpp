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

Flow::Flow(double A, double B, double C) {
	type_ = ABC;
	A_ = A;
	B_ = B;
	C_ = C;
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
		case ABC:
			switch(dim) {
				case 0:
					v = A_*sin(z) + C_*cos(y);
					break;
				case 1:
					v = B_*sin(x) + A_*cos(z);
					break;
				case 2:
					v = C_*sin(y)+B_*cos(x);
					break;
				default:
					ErrorMsg("Flow::velocity: dimension error");
			}
			break;
		default:
			ErrorMsg("Flow::velocity: Flow type not found");
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
		case ABC:
			switch(dim) {
				case 1: // \partial_y u_x
					g = -C_*sin(y);
					break;
				case 2: // \partial_z u_x
					g = A_*cos(z);
					break;
				case 3: // \partial_x u_y
					g = B_*cos(x);
					break;
				case 5: // \partial_z u_y
					g = -A_*sin(z);
					break;
				case 6: // \partial_x u_z
					g = -B_*sin(x);
					break;
				case 7: // \partial_y u_z
					g = C_*cos(y);
					break;
				default: // \partial_x u_x, \partial_y u_y, \partial_z u_z
					g = 0;
			}
			break;
		default:
			ErrorMsg("Flow::gradient: Flow type not found");
			break;
	}
	return g;
}
