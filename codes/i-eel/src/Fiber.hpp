//
//  Fiber.hpp
//  i-eel
//
//  Created by Jeremie Bec on 25/11/2020.
//

#ifndef Fiber_hpp
#define Fiber_hpp

#include "basics/Utilities.h"
#include "Flow.hpp"

using namespace std;
using namespace arma;

typedef Mat<double> MyMat;
typedef Col<double> MyCol;
typedef SpMat<double> MySpMat;

class Fiber {
public:
	// Calss initialisation
	Fiber();
	Fiber(int Ns, double L, double mu, double beta);
	Fiber(int N, double L, double mu, double beta, Flow& U, std::default_random_engine& rng);
	Fiber(int N, double L, double mu, double beta, Flow& U, vector<double> p);
	void read(int N, double L, double mu, double beta, Flow& U, const std::string& filebase);
	void alloc();
	// save the fiber configuration to a file
	void save(const std::string& filebase) const;
	// set the forcing parameters
	void setforcing(double k, double omega, double A);
	// compute the tension forces
	void calc_tension();
	// time evolution over a time step dt
	void evol(double dt, Flow& U);
	// compute the end-to-end length
	inline double endtoend() const {
		return sqrt((X_(Ns_,0)-X_(0,0))*(X_(Ns_,0)-X_(0,0)) + (X_(Ns_,1)-X_(0,1))*(X_(Ns_,1)-X_(0,1)) + (X_(Ns_,2)-X_(0,2))*(X_(Ns_,2)-X_(0,2)));
	};
	// compute the average stretching (in principle equal to 1
	inline double meanstretch() const {
		return sum( D1X_.col(0)%D1X_.col(0) + D1X_.col(1)%D1X_.col(1) + D1X_.col(2)%D1X_.col(2) )/(double)Ns_;
	};
	// compute the maximal curvature
	inline double maxcurv() const {
		MyCol K = D2X_.col(0)%D2X_.col(0) + D2X_.col(1)%D2X_.col(1) + D2X_.col(2)%D2X_.col(2);
		return sqrt(K.max());
	};
	// save the position and the tension in ascii
	void saveascii(const std::string& filebase) {
		X_.save(filebase+"_X.dat",raw_ascii);
		T_.save(filebase+"_T.dat",raw_ascii);
	}
	// set the forcing parameters
	
	// return internal variables
	inline double L() const { return L_; };
	inline int N() const { return Ns_; };
	
private:
	// Physical parameters
	double L_;
	double mu_;
	// Forcing parameters (set by default to 0)
	double k_ = 0.0;
	double om_ = 0.0;
	double A_ = 0.0;
	// Numerical parameters
	int Ns_;
	double ds_;
	double beta_;
	double dt_;
	// Internal variables
	MyMat X_, Xold_, Gold_; // position
	MyMat D1X_, D2X_, D3X_, D4X_; // derivatives
	MyCol NormXi_, T_, D1T_; // tension
	MyMat  Uf_, D1Uf_; // fluid velocity at particle position
	MyMat F_; // forcing
	// Derivative matrices
	MySpMat D1_, D2_, D3_, D4_, LapDirich_, Op4_;
	// Internal functions
	void diff();
	void interp_U(Flow& U);
	void calc_force();
	inline double wrap(double x, double L) const { return x-L*floor(x/L); };
};

#endif /* Fiber_hpp */
