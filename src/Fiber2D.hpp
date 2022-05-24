//
//  Fiber2D.hpp
//  i-eel
//
//  Created by Jeremie Bec on 19/01/2021.
//

#ifndef Fiber2D_hpp
#define Fiber2D_hpp

#include "basics/Utilities.h"
#include "Flow2D.hpp"

using namespace std;
using namespace arma;

typedef Mat<double> MyMat;
typedef Col<double> MyCol;
typedef SpMat<double> MySpMat;

MyMat readQ(const std::string& filebase);
void set_seed(void);
class Fiber2D {
public:
    // Calss initialisation
    Fiber2D();
    Fiber2D(int Ns, double L, double zeta, double E, double beta);
    Fiber2D(int N, double L, double zeta, double E, double beta, Flow2D& U, std::default_random_engine& rng);
    Fiber2D(int N, double L, double zeta, double E, double beta, Flow2D& U, vector<double> p);
    void read(Flow2D& U, const std::string& filebase, double t = 0, double dt = 0);
    void read(int N, double L, double zeta, double E, double beta, Flow2D& U, const std::string& filebase, vector<double>, double, double);
    void set_seed(void);
    void alloc();
    // save the fiber configuration to a file
    void save(const std::string& filebase, Flow2D&U);
    // set the forcing parameters
    void setforcing(int k, double om);
    void setforcing(vector<double> p, double A);
    // compute the tension forces
    void calc_tension();
    // time evolution over a time step dt
    void evol(double dt, Flow2D& U);
    // compute the position of the center of mass
    inline double getcenter(int dim) {
        double x = 0;
        for(int is=0; is<=Ns_; ++is)
        x += X_(is,dim);
        return x/((double)(Ns_+1));
    };
    // compute the mean translational velocity
    inline double getvelocity(int dim) {
        double v = 0;
        for(int is=0; is<=Ns_; ++is)
        v += X_(is,dim)-Xold_(is,dim);
        return v/((double)(Ns_+1))/dt_;
    };
    // compute the end-to-end length
    inline double endtoend() const {
        return sqrt((X_(Ns_,0)-X_(0,0))*(X_(Ns_,0)-X_(0,0)) + (X_(Ns_,1)-X_(0,1))*(X_(Ns_,1)-X_(0,1)) );
    };
    // compute the average stretching (in principle equal to 1)
    inline double meanstretch() const {
        return sum( D1X_.col(0)%D1X_.col(0) + D1X_.col(1)%D1X_.col(1) )/(double)Ns_;
    };
    // compute the maximal curvature
    inline double maxcurv() const {
        MyCol K = D2X_.col(0)%D2X_.col(0) + D2X_.col(1)%D2X_.col(1) ;
        return sqrt(K.max());
    };
    // save the position and the tension in ascii
    void saveascii(const std::string& filebase) {
        X_.save(filebase+"_X.dat",raw_ascii);
        T_.save(filebase+"_T.dat",raw_ascii);
    };
    double wind(Flow2D& U);
    vector<double> orientation();
    int calc_buckle();
    
    // return internal variables
    inline double L() const { return L_; };
    inline int N() const { return Ns_; };
    
private:
    // Time
    double t_ = 0.0;
    // Physical parameters
    double L_ = 1.0 ;    // Fiber length
    double zeta_; // Viscous friction coefficient
    double E_;    // Young modulus
    // Forcing parameters (set by default to 0)
    vector<double> p_; // direction of the helix
    double nu_ = 0.0;  // wavenumber
    double om_ = 0.0;  // frequency
    double A_ = 0.0;   // internal amplitude
    // Numerical parameters
    int Ns_;
    double ds_;
    double beta_;
    double dt_;
    vector<double> Rtot_, Rew_, State_, Orient_, Vent_;
    vector<int> Boucle_;
    // Internal variables
    MyMat X_, Xold_, Gold_, X_Qlearning_; // position
    MyMat D1X_, D2X_, D3X_, D4X_; // derivatives
    MyCol NormXi_, T_, D1T_; // tension
    MyMat  Uf_, D1Uf_; // fluid velocity at particle position
    MyMat F_, D1F_; // forcing
    
    // Derivative matrices
    MySpMat D1_, D2_, D3_, D4_, LapDirich_, Op4_;
    // Internal functions
    void diff();
    void interp_U(Flow2D& U);
    void calc_force();
    inline double wrap(double x, double L) const { return x-L*floor(x/L); };
};

#endif /* Fiber2D_hpp */
