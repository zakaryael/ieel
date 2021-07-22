//
//  learn2D_cellflow.cpp
//  i-eel
//
//  Created by Jeremie Bec on 27/05/2021.
//

#include "src/Fiber2D.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"

using namespace std;

int main(int argc, char* argv[]) {
    WriteProcessInfo(argc, argv);
    
    string purpose("Q-Learning of a 2D fiber with no flow.");
    
    ArgList args(argc, argv, purpose);
    
    // Read the command-line options
    args.section("Program options");
    const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
    const string indir = args.getpath("-in", "--indir", "data/", "input directory");
    const double L = args.getreal("-L", "--length", 1.0, "fiber length");
    const double zeta = args.getreal("-z", "--zeta", 1e5, "friction coefficient");
    const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
    const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
    const double Tmax = args.getreal("-T", "--time",10000.0, "integration time");
    const int Nout = args.getint("-nout", "--step_out", 2000, "output period (in number of timesteps)");
    const double dt = args.getreal("-dt", "--timestep", 0.0005, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
    const double om = args.getreal("-om", "--frequency", 2, "Forcing frequency");
    const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
    const int Nlearning = args.getint("-nl", "--step_learning", 200, "learning period (in number of timesteps)");
    const double gamma = args.getreal("-gamma", "--discountrate", 0.9995, "Discount rate");
    const double learnrate = args.getreal("-lr", "--learnrate", 0.005, "Learning rate");
    const double epsil = args.getreal("-eps", "--epsilon", 0.0, "Rate of random exploration");
    const double qinit = args.getreal("-q0", "--qinit", 0.25, "Initial Q entries");
    const string inQ = args.getstr("-Qinit", "--initialQ", "", "input file from which Q is read");
    
    args.check();
    mkdir(outdir);
    args.save(outdir);
    
    // fluid flow set to 0
    Flow2D U(Null);
    // define the fiber
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"Generating a straight fiber of length "<<L<<endl;
    cout<<"with fric. coeff. "<<zeta<<" and Young modulus"<<E<<endl;
    vector<double> p(2);
    p.at(0) = -1.0;
    p.at(1) = 0.0;
    
    cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
    Fiber2D Fib(Ns,L,zeta,E,beta,U,p);
    
    // set the forcing
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
    Fib.setforcing(k, om);
    MyMat Q;
    Q.set_size(12,8);
    if(strcmp(inQ.c_str(),"")==0) {
        cout<<"Start Q from scratch"<<endl;
        for (int i = 0; i <12; i++) {
            for (int j = 0; j<8; j++)
            Q(i,j) = qinit;
        }
    }
    else {
        cout<<"Restart Q from file "<<inQ<<endl;
        Q = readQ(inQ);
    }
    MyCol Ampl;
    double a0 = zeta*alpha*om/(2.0*M_PI*(double)k/L);
    Ampl.set_size(4);
    Ampl(0) = 0;
    Ampl(1) = a0/3.0;
    Ampl(2) = 2.0*a0/3.0;
    Ampl(3) = a0;

    Fib.initQlearning(Q,gamma,learnrate,1.0,Ampl,U,epsil);
    
    // time loop
    cout<<endl<<"------------------------------------------------"<<endl;
    int nstep = round(Tmax/dt);
    cout<<"Starting the time loop over "<<nstep<<" steps"<<endl;
    cout<<"output every "<<Nout<<" steps"<<endl;
    int it = 0;
    double t = 0;
    
    //generator.seed((unsigned) time(&tt));
    //set_seed();
    while(it<nstep) {
        //generator.seed((unsigned) time(&tt));
        if((it % Nout)==0) {
            cout<<setprecision(4);
            cout << showpoint;
            cout<<"t = "<<t<<setw(10)<<"X = "<<Fib.getcenter(0)<<setw(10)<<"Vx = "<<Fib.printreward()<<" Action: "<< Fib.getaction()<<" State: "<< Fib.getstate()<<endl;
            Fib.save(outdir+"fiber"+i2s(it)+".nc",U);
        }
        Fib.evol(dt,U);
        t += dt;
        it++;
        if((it % Nlearning)==0) {
            Fib.setepsilon(Fib.epsilon()/(double)it);
            Fib.Qupdate(U);
        }
    }
    
    return 1;
}

