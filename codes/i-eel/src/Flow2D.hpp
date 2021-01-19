//
//  Flow2D.hpp
//  i-eel
//
//  Created by Jeremie Bec on 19/01/2021.
//

#ifndef Flow2D_hpp
#define Flow2D_hpp

#include "basics/Utilities.h"

enum flow2Dtype { Null, Shear, ABC };

class Flow2D {
public:
    Flow2D() {type_ = Null; };
    Flow2D(double sig);
    double velocity(double x, double y, int dim);
    double gradient(double x, double y, int dim);
private:
    flow2Dtype type_ = Null;
    double *shear_; // contains the mean gradient \partial_j u_i in element 2*i+j
};


#endif /* Flow2D_hpp */
