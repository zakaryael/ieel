//
//  test_buckling2D.cpp
//  i-eel
//
//  Created by Jeremie Bec on 19/01/2021.
//

#include "src/Fiber2D.hpp"
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
    const double sig = args.getreal("-sig", "--shear", 1.0, "shear amplitude");
    const double dev = args.getreal("-d", "--deviation", 0.1, "initial deviation of the orientation");
    const double L = args.getreal("-L", "--length", 1.0, "fiber length");
    const double zeta = args.getreal("-z", "--zeta", 1e6, "friction coefficients");
    const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
    const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
    const double Tmax = args.getreal("-T", "--time", 25.0, "integration time");
    const int Nout = args.getint("-nout", "--step_out", 200, "output period (in number of timesteps)");
    const double dt = args.getreal("-dt", "--timestep", 1e-3, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    args.check();
    mkdir(outdir);
    args.save(outdir);
    
    // define the flow
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"Generating a pure shear flow"<<endl;
    Flow2D U(sig);
    cout<<"     [\t\t"<<U.gradient(0,0,0)<<"\t"<<U.gradient(0,0,1)<<"\t]"<<endl;
    cout<<"DU = [\t\t"<<U.gradient(0,0,2)<<"\t"<<U.gradient(0,0,3)<<"\t]"<<endl;
    
    // define the fiber
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"Generating a straight fiber of length "<<L<<endl;
    cout<<"with fric. coeff. "<<zeta<<" and Young modulus"<<E<<endl;
    vector<double> p(3);
    p.at(0) = sqrt(1-dev*dev);
    p.at(1) = -dev;
    cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
    Fiber2D Fib(Ns,L,zeta,E,beta,U,p);
    
    // time loop
    cout<<endl<<"------------------------------------------------"<<endl;
    int nstep = round(Tmax/dt);
    cout<<"Starting the time loop over "<<nstep<<" steps"<<endl;
    cout<<"output every "<<Nout<<" steps"<<endl;
    int it = 0;
    double t = 0;
    
    while(it<nstep) {
        if((it % Nout)==0) {
            cout<<setprecision(4);
            cout << showpoint;
            cout<<"t = "<<t<<setw(10)<<"Lee = "<<Fib.endtoend()<<endl;
            Fib.save(outdir+"fiber"+i2s(it)+".ff",U);
        }
        Fib.evol(dt,U);
        t += dt;
        it++;
     }
    Fib.save(outdir+"fiber"+i2s(it)+".ff",U);
    return 1;
}
