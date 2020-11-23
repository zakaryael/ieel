//
//  fibers.h
//  i-eel
//
//  Created by Jeremie Bec on 23/11/2020.
//

#ifndef fibers_h
#define fibers_h

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/basisfunc.h"
#include "channelflow/cfmpi.h"
#include "channelflow/chebyshev.h"
#include "channelflow/realprofile.h"
#include "channelflow/realprofileng.h"
#include "channelflow/flowfield.h"
#include "channelflow/diffops.h"

#include <random>
#include <armadillo>

using namespace std;
using namespace chflow;
using namespace arma;

typedef Mat<Real> MyMat;
typedef Col<Real> MyCol;
typedef SpMat<Real> MySpMat;

class Fiber {
public:
	Fiber();
	Fiber(int Ns, Real L, Real mu, Real beta, Real gamma, CfMPI* cfmpi=NULL);
	void init(int N, Real L, Real mu, Real beta, Real gamma, FlowField& u, FlowField& du, std::default_random_engine& rng,  CfMPI* cfmpi=NULL);
	void read(int N, Real L, Real mu, Real beta, Real gamma, const std::string& filebase, int ifib, FlowField& u, FlowField& du, CfMPI* cfmpi=NULL);
	void alloc();
	void save(const std::string& filebase) const;
	void calc_tension();
	void evol(Real dt, FlowField& u, FlowField& du);
	
	inline Real endtoend() const {
		return sqrt((X_(Ns_,0)-X_(0,0))*(X_(Ns_,0)-X_(0,0)) + (X_(Ns_,1)-X_(0,1))*(X_(Ns_,1)-X_(0,1)) + (X_(Ns_,2)-X_(0,2))*(X_(Ns_,2)-X_(0,2)));
	};
	inline Real meanstretch() const {
		return sum( D1X_.col(0)%D1X_.col(0) + D1X_.col(1)%D1X_.col(1) + D1X_.col(2)%D1X_.col(2) )/(Real)Ns_;
	};
	inline Real maxcurv() const {
		MyCol K = D2X_.col(0)%D2X_.col(0) + D2X_.col(1)%D2X_.col(1) + D2X_.col(2)%D2X_.col(2);
		return sqrt(K.max());
	};
	void saveascii(const std::string& filebase) {
		X_.save(filebase+"_X.dat",raw_ascii);
		T_.save(filebase+"_T.dat",raw_ascii);
	}
	inline Real L() const { return L_; };
	inline int N() const { return Ns_; };
	
private:
	Real L_;
	Real mu_;
	int Ns_;
	Real ds_;
	Real beta_;
	Real gamma_;
	
	Real dt_;
	
	MyMat X_, Xold_, Gold_;
	MyMat D1X_, D2X_, D3X_, D4X_;
	MyCol NormXi_, T_, D1T_;
	MyMat  Uf_, D1Uf_;
	
	MySpMat D1_, D2_, D3_, D4_, LapDirich_, Op4_;
	
	CfMPI* cfmpi_ = nullptr;
	
	void diff();
	void interp_U(FlowField& u, FlowField& du);
	inline Real wrap(Real x, Real L) const { return x-L*floor(x/L); };
};

#endif /* fibers_h */
