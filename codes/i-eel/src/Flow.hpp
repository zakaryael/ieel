//
//  Flow.hpp
//  i-eel
//
//  Created by Jeremie Bec on 25/11/2020.
//

#ifndef Flow_hpp
#define Flow_hpp

#include "basics/Utilities.h"

enum flowtype { Null, Shear, ABC };

class Flow {
public:
	Flow() {type_ = Null; };
	Flow(double sig);
	double velocity(double x, double y, double z, int dim);
	double gradient(double x, double y, double z, int dim);
private:
	flowtype type_ = Null;
	double *shear_; // contains the mean gradient \partial_j u_i in element 3*i+j
};

#endif /* Flow_hpp */
