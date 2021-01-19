//
//  Flow2D.cpp
//  i-eel
//
//  Created by Jeremie Bec on 19/01/2021.
//

#include "Flow2D.hpp"

Flow2D::Flow2D(double sig) {
    type_ = Shear;
    shear_ = new double[4];
    for(int idim=0; idim<4; ++idim)
    shear_[idim] = 0;
    shear_[1] = sig; // \partial_y u_x = sig
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
        default:
            ErrorMsg("Flow type not found");
            break;
    }
    return g;
}
