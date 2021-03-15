//
//  swim2D_learned_cellflow.cpp
//  i-eel
//
//  Created by Jeremie Bec on 09/03/2021.
//
#include "src/Fiber2D.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"

using namespace std;

int main(int argc, char* argv[]) {
	WriteProcessInfo(argc, argv);
	
	string purpose("Naive Swimming of a 2D fiber in a cellular flow.");
	
	ArgList args(argc, argv, purpose);
	
	// Read the command-line options
	args.section("Program options");
	const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
	const string indir = args.getpath("-in", "--indir", "data/", "input directory");
	const double L = args.getreal("-L", "--length", 1.0, "fiber length");
	const double zeta = args.getreal("-z", "--zeta", 1e5, "friction coefficient");
	const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
	const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
	const double Tmax = args.getreal("-T", "--time",600.0, "integration time");
	const int Nout = args.getint("-nout", "--step_out", 400, "output period (in number of timesteps)");
	const double dt = args.getreal("-dt", "--timestep", 0.0005, "time step");
	const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
	const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
	const double om = args.getreal("-om", "--frequency", 2, "Forcing frequency");
	const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
	const double u = args.getreal("-U", "--Velocity", 0.05, "Velocity amplitude");
	const int Nlearning = args.getint("-nl", "--step_learning", 1, "learning period (in number of timesteps)");
	const double u0 = args.getreal("-slim", "--speed", 0.01, "Vitesse limite");
	
	args.check();
	mkdir(outdir);
	args.save(outdir);
	
	// fluid flow set to 0
	Flow2D U(Cellular);
	U.initcellular(u);
	// define the fiber
	cout<<endl<<"------------------------------------------------"<<endl;
	cout<<"Generating a straight fiber of length "<<L<<endl;
	cout<<"with fric. coeff. "<<zeta<<" and Young modulus"<<E<<endl;
	vector<double> p(2);
	time_t tt;
	std::default_random_engine rng;
	rng.seed((unsigned) time(&tt));
	std::uniform_real_distribution<double> Unif(0.0,1.0);
	p.at(0) = Unif(rng);
	p.at(1) = Unif(rng);
	double r2;
	while((r2 = p.at(0)*p.at(0)+p.at(1)*p.at(1))>1) {
		p.at(0) = Unif(rng);
		p.at(1) = Unif(rng);
	}
	p.at(0) /= sqrt(r2);
	p.at(1) /= sqrt(r2);
	if(p.at(0)>0)
		p.at(0) *= -1;
	
	cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
	Fiber2D Fib(Ns,L,zeta,E,beta,U,p);
	
	// set the forcing
	cout<<endl<<"------------------------------------------------"<<endl;
	cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
	Fib.setforcing(k, om);
	MyMat Q;
	Q.set_size(12,8);
	for(int i=0; i<12; ++i) {
		for(int j=0; j<4; ++j)
		Q(i,j) = 0;
	}
	// learned strategy from the best performing learning over an ensemble of 50 trajectories
	for(int i=0; i<6; ++i) // States corresponding to being wrongly oriented
	Q(i,2) = 10;
	for(int i=6; i<12; ++i) // States corresponding to being rightly oriented
	Q(i,3) = 10;
	MyCol Ampl;
	double a0 = zeta*alpha*om/(2.0*M_PI*(double)k/L);
	Ampl.set_size(2);
	Ampl(0) = 0;
	Ampl(1) = a0;
	
	Fib.initQlearning(Q,0,0,u0,Ampl,U);
	
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
			cout<<"t = "<<t<<setw(10)<<"Vx = "<<Fib.getvelocity(0)<< " Action : " << Fib.getaction() << endl;
			x0 = x;
			Fib.save(outdir+"fiber"+i2s(it)+".nc",U);
		}
		Fib.evol(dt,U);
		t += dt;
		it++;
		if((it % Nlearning)==0)
			Fib.Qupdate(U);
	}
	
	return 1;
}
