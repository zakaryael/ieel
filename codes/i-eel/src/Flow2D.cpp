//
//  Flow2D.cpp
//  i-eel
//
//  Created by Jeremie Bec on 19/01/2021.
//

#include "Flow2D.hpp"

Flow2D::Flow2D(double sig) {
	initshear(sig);
}

Flow2D::Flow2D(flow2Dtype type) {
	if(type==Shear)
		Flow2D::initshear(1.0);
	else
		type_ = type;
}

void Flow2D::initshear(double sig) {
	type_ = Shear;
	shear_ = new double[4];
	for(int idim=0; idim<4; ++idim)
	shear_[idim] = 0;
	shear_[1] = sig; // \partial_y u_x = sig
}

void Flow2D::initcellular(double u) {
	type_ = Cellular;
	u_ = u;
}

double Flow2D::velocity(double x, double y, int dim) {
    double v = 0.0;
    assert(dim>=0 && dim<2);
    switch (type_) {
        case Null:
            break;
        case Shear:
            v = shear_[2*dim]*x+shear_[2*dim+1]*y;
            break;
		case Cellular:
			// Stream function: cos(2\pi x) cos(2\pi y) / (2\pi)
			switch(dim) {
				case 0:
					v = -u_*cos(M_PI*x/L_)*sin(M_PI*y/L_);
					break;
				case 1:
					v = u_*sin(M_PI*x/L_)*cos(M_PI*y/L_);
					break;
				default:
					ErrorMsg("Flow2D::velocity: error trying to access an unexisting dimension");
					break;
			}
			break;
        default:
            ErrorMsg("Flow type not found");
            break;
    }
    return v;
}

double Flow2D::gradient(double x, double y, int dim) {
    double g = 0.0;
    assert(dim>=0 && dim<4);
    switch (type_) {
        case Null:
            break;
        case Shear:
            g = shear_[dim];
            break;
		case Cellular:
			switch (dim) {
				case 0: // \partial_x u_x
					g = (M_PI/L_)*u_*sin(M_PI*x/L_)*sin(M_PI*y/L_);
					break;
				case 1: // \partial_y u_x
					g = -(M_PI/L_)*u_*cos(M_PI*x/L_)*cos(M_PI*y/L_);
					break;
				case 2: // \partial_x u_y
					g = (M_PI/L_)*u_*cos(M_PI*x/L_)*cos(M_PI*y/L_);
					break;
				case 3: // \partial_y u_y
					g = -(M_PI/L_)*u_*sin(M_PI*x/L_)*sin(M_PI*y/L_);
					break;
				default:
					ErrorMsg("Flow2D::gradient: error trying to access an unexisting dimension");
			}
			break;
        default:
            ErrorMsg("Flow type not found");
            break;
    }
    return g;
}
