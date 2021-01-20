//
//  test_swim2D.cpp
//  i-eel
//
//  Created by Jeremie Bec on 20/01/2021.
//

#include "src/Fiber2D.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"

using namespace std;

int main(int argc, char* argv[]) {
	WriteProcessInfo(argc, argv);
	
	string purpose("Test on the swimming of a 2D fiber in the absence of flow.");
	
	ArgList args(argc, argv, purpose);
	
	// Read the command-line options
	args.section("Program options");
	const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
	const double L = args.getreal("-L", "--length", 1.0, "fiber length");
	const double zeta = args.getreal("-z", "--zeta", 1e5, "friction coefficient");
	const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
	const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
	const double Tmax = args.getreal("-T", "--time", 50.0, "integration time");
	const int Nout = args.getint("-nout", "--step_out", 200, "output period (in number of timesteps)");
	const double dt = args.getreal("-dt", "--timestep", 1e-3, "time step");
	const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
	const int k = args.getint("-k", "--wavenumber", 5, "Forcing wavenumber");
	const double om = args.getreal("-om", "--frequency", 5.0, "Forcing frequency");
	const double alpha = args.getreal("-alpha", "--alpha", 0.5, "Force amplitude");
	
	args.check();
	mkdir(outdir);
	args.save(outdir);
	
	// fluid flow set to 0
	Flow2D U(0);
	
	// define the fiber
	cout<<endl<<"------------------------------------------------"<<endl;
	cout<<"Generating a straight fiber of length "<<L<<endl;
	cout<<"with fric. coeff. "<<zeta<<" and Young modulus"<<E<<endl;
	vector<double> p(2);
	p.at(0) = 1;
	p.at(1) = 0;
	cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
	Fiber2D Fib(Ns,L,zeta,E,beta,U,p);
	
	// set the forcing
	cout<<endl<<"------------------------------------------------"<<endl;
	cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
	cout<<"in the direction ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
	Fib.setforcing(p, k, om, alpha);
	
	// time loop
	cout<<endl<<"------------------------------------------------"<<endl;
	int nstep = round(Tmax/dt);
	cout<<"Starting the time loop over "<<nstep<<" steps"<<endl;
	cout<<"output every "<<Nout<<" steps"<<endl;
	int it = 0;
	double t = 0;
	
	double x0,x;
	x0 = Fib.getcenter(0);
	while(it<nstep) {
		if((it % Nout)==0) {
			cout<<setprecision(4);
			cout << showpoint;
			x = Fib.getcenter(0);
			cout<<"t = "<<t<<setw(10)<<"Vx = "<<(x-x0)/(Nout*dt)<<endl;
			x0 = x;
			Fib.save(outdir+"fiber"+i2s(it)+".nc",U);
		}
		
		Fib.evol(dt,U);
		t += dt;
		it++;
	}
	Fib.save(outdir+"fiber"+i2s(it)+".nc",U);
	return 1;
}
