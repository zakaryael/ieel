//
//  Flow2D.hpp
//  i-eel
//
//  Created by Jeremie Bec on 19/01/2021.
//

#ifndef Flow2D_hpp
#define Flow2D_hpp

#include "basics/Utilities.h"

enum flow2Dtype { Null, Shear, Cellular };

class Flow2D {
public:
    Flow2D() {type_ = Null; };
	Flow2D(flow2Dtype type);
	Flow2D(double sig);
	void initshear(double sig);
	void initcellular(double u, double L = 1.0);
	double velocity(double x, double y, int dim);
    double gradient(double x, double y, int dim);
private:
    flow2Dtype type_ = Null;
    double *shear_; // contains the mean gradient \partial_j u_i in element 2*i+j
	double u_;
	double L_;
};


#endif /* Flow2D_hpp */
