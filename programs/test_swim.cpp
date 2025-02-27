//
//  test_swim.cpp
//  i-eel
//
//  Created by Jeremie Bec on 28/11/2020.
//

#include "src/Fiber.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"

using namespace std;

int main(int argc, char* argv[]) {
	WriteProcessInfo(argc, argv);
	
	string purpose("Test on the development of the buckling instability in a pure shear flow.");
	
	ArgList args(argc, argv, purpose);
	
	// Read the command-line options
	args.section("Program options");
	const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
	const double L = args.getreal("-L", "--length", 2.0*M_PI, "fiber length");
	const double zeta = args.getreal("-z", "--zeta", 1, "friction coefficient");
	const double E = args.getreal("-E", "--EI", 1e-3, "Young modulus");
	const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
	const double Tmax = args.getreal("-T", "--time", 50.0, "integration time");
	const int Nout = args.getint("-nout", "--step_out", 200, "output period (in number of timesteps)");
	const double dt = args.getreal("-dt", "--timestep", 1e-3, "time step");
	const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
	const int k = args.getint("-k", "--wavenumber", 5, "Forcing wavenumber");
	const double om = args.getreal("-om", "--frequency", 5.0, "Forcing frequency");
	//const double Tw = args.getreal("-tw", "--twist", 0.1, "Forcing twist");
	const int eps = args.getint("-e", "--chiral", 1, "Forcing chirality");
    const double alpha = args.getreal("-alpha", "--alpha", 0.5, "Force amplitude");
	
	args.check();
	mkdir(outdir);
	args.save(outdir);
	
	// fluid flow set to 0
	Flow U(1,1,1);
	
	// define the fiber
	cout<<endl<<"------------------------------------------------"<<endl;
	cout<<"Generating a straight fiber of length "<<L<<endl;
	cout<<"with fric. coeff. "<<zeta<<" and Young modulus"<<E<<endl;
	vector<double> p(3);
	p.at(0) = 1;
	p.at(1) = 0;
	p.at(2) = 0;
	cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<","<<p.at(2)<<")"<<endl;
	Fiber Fib(Ns,L,zeta,E,beta,U,p);
	
	// set the forcing
	p.at(0) = 1/sqrt(3);
	p.at(1) = 1/sqrt(3);
	p.at(2) = 1/sqrt(3);
	cout<<endl<<"------------------------------------------------"<<endl;
	cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
	cout<<"in the direction ("<<p.at(0)<<","<<p.at(1)<<","<<p.at(2)<<")"<<endl;
	double nu = 2.0*M_PI*(double)k/L;
    double Rdeb = (1.0/nu)*((1.0/alpha) - sqrt((1.0/(alpha*alpha))-1.0));
	cout<<"Predicted velocity: "<<-(alpha*om*Rdeb*sqrt(1.0-nu*nu*Rdeb*Rdeb)/2.0)<<endl;
	Fib.setforcing2(p, nu, om, alpha, eps);
	
	// time loop
	cout<<endl<<"------------------------------------------------"<<endl;
	int nstep = round(Tmax/dt);
	cout<<"Starting the time loop over "<<nstep<<" steps"<<endl;
	cout<<"output every "<<Nout<<" steps"<<endl;
	int it = 0;
	double t = 0;
	
	double x0, x, y0, y, z, z0;
	x0 = Fib.getcenter(0);
	y0 = Fib.getcenter(1);
	z0 = Fib.getcenter(2);
	while(it<nstep) {
		if((it % Nout)==0) {
			cout<<setprecision(4);
			cout << showpoint;
			x = Fib.getcenter(0);
			y = Fib.getcenter(1);
			z = Fib.getcenter(2);
			cout<<"t = "<<t<<setw(10)<<"V = "<<-(p.at(0)*(x-x0)+p.at(1)*(y-y0)+p.at(2)*(z-z0))/(Nout*dt)<<endl;
			x0 = x;
			y0 = y;
			z0 = z;
			Fib.save(outdir+"fiber"+i2s(it)+".nc",U);
		}
		
		Fib.evol(dt,U);
		t += dt;
		it++;
	}
	Fib.save(outdir+"fiber"+i2s(it)+".nc",U);
	return 1;
}
