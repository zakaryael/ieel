//
//  learn2D_reduced.cpp
//  i-eel
//
//  Created by Jeremie Bec on 30/07/2021.
//

#include "src/Fiber2D.hpp"
#include "src/QLearning.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"

using namespace std;

int main(int argc, char* argv[]) {
    WriteProcessInfo(argc, argv);
    
    string purpose("Q-Learning of a 2D fiber in a cellular flow.");
    
    ArgList args(argc, argv, purpose);
    
    // Read the command-line options
    args.section("Program options");
    const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
    const string indir = args.getpath("-in", "--indir", "", "input directory to a binary file from which initial Q is read");
    const string Qinitdir = args.getpath("-Qdir", "--Qinitdir", "", "input directory to a csv file from which the initial Q is read");
    const string Pinitdir = args.getpath("-Pidir", "--Pinitdir", "", "input directory to a csv file from which the initial policy is read");
    const double L = args.getreal("-L", "--length", 1.0, "fiber length");
    const double zeta = args.getreal("-z", "--zeta", 5e4, "friction coefficient");
    const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
    const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
    const double Tmax = args.getreal("-T", "--time",20000.0, "integration time");
    const int Nout = args.getint("-nout", "--step_out", 2000, "output period (in number of timesteps)");
    const double dt = args.getreal("-dt", "--timestep", 0.0005, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
    const double om = args.getreal("-om", "--frequency", 2, "Forcing frequency");
    const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
    const double u = args.getreal("-U", "--Velocity", 0.05, "Velocity amplitude");
    const int Nlearning = args.getint("-nl", "--step_learning", 200, "learning period (in number of timesteps)");
    const int Noutlearning = args.getint("-nlout", "--step_save_learning", 200, "learning output period (in number of timesteps)");
    const double gamma = args.getreal("-gamma", "--discountrate", 0.9995, "Discount rate");
    const double learnrate = args.getreal("-lr", "--learnrate", 0.005, "Learning rate");
    const double epsil = args.getreal("-eps", "--epsilon", 0.0, "Rate of random exploration");
    const double u0 = args.getreal("-slim", "--speed", 0.01, "Vitesse limite");
    const double qinit = args.getreal("-q0", "--qinit", 2, "Initial Q entries");
    const int learning = args.getint("-lrn", "--learning", 1, "Input 1 for swimming with learning 0 otherwise");
    const int noflow = args.getint("-nfl", "--noflow", 0, "Input 1 for no flow 0 for cellular flow");
    args.check();
    mkdir(outdir);
    args.save(outdir);
    
    // fluid flow set to cellular with amplitude u
    Flow2D U(Cellular);
    U.initcellular(u);
    // To set the fluid flow set to 0:
    if (noflow == 1){Flow2D U(Null);}
    
    // define the fiber
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"Generating a straight fiber of length "<<L<<endl;
    cout<<"with fric. coeff. "<<zeta<<" and Young modulus"<<E<<endl;
    vector<double> p(2);
    time_t tt;
    std::default_random_engine rng;
    rng.seed((unsigned) time(&tt));
    std::uniform_real_distribution<double> Unif(0.0,1.0);
    /*p.at(0) = Unif(rng);
    p.at(1) = Unif(rng);
    double r2;
    while((r2 = p.at(0)*p.at(0)+p.at(1)*p.at(1))>1) {
        p.at(0) = Unif(rng);
        p.at(1) = Unif(rng);
    }
    p.at(0) /= sqrt(r2);
    p.at(1) /= sqrt(r2);
    if(p.at(0)>0)
        p.at(0) *= -1; */
    p.at(0) = -1;
    p.at(1) = 0;
    cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
    Fiber2D Fib(Ns,L,zeta,E,beta,U,p);
    // set the forcing
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
    Fib.setforcing(k, om);
    
    // define the learning
    MyMat Q;
    MyMat Pi;
    Pi.set_size(4, 3);
    //Pi.ones();
    Pi(0,0) = 1; Pi(0,1) = 0; Pi(0,2) = 0;
    Pi(1,0) = 1; Pi(1,1) = 0; Pi(1,2) = 0;
    Pi(2,0) = 1; Pi(2,1) = 0; Pi(2,2) = 0;
    Pi(3,0) = 0; Pi(3,1) = 0; Pi(3,2) = 1;
    
//    if(strcmp(Qinitdir.c_str(),"")==0){
//        Q.load("Q.csv", arma::csv_ascii);
//    }
//    else{
//        Q.load(Qinitdir, arma::csv_ascii);
//    }
    if(strcmp(indir.c_str(),"")==0) {
        cout<<"Start Q from scratch..."<<qinit<<endl;
        Q.set_size(4,3);
         for (int i = 0; i<4; i++) {
             for (int j = 0; j<3; j++)
                 Q(i,j) = qinit;
         }
    }
    else {
        cout<<"Start Q from file "<<indir<<"learn.bin"<<endl;
        Q = readlastQ(indir+"learn.bin",4,3);
    }
    
    MyCol Ampl;
    double a0 = zeta*alpha*om/(2.0*M_PI*(double)k/L);
    Ampl.set_size(2);
    Ampl(0) = 0;
    Ampl(1) = a0;
    QLearning QL(Q, Pi, gamma, learnrate, u0, Ampl, epsil, true);
    char cname[512];
    string fname = outdir+"learn.bin";
    strcpy(cname, fname.c_str());
    FILE *fout = fopen(cname,"w");
    fclose(fout);
    
    // time loop
    cout<<endl<<"------------------------------------------------"<<endl;
    int nstep = round(Tmax/dt);
    cout<<"Starting the time loop over "<<nstep<<" steps"<<endl;
    cout<<"output every "<<Nout<<" steps"<<endl;
    int it = 0;
    double t = 0;
    
    
    int state, previous_state = QL.compute_state(Fib.wind(U), Fib.orientation()); //computing s0
    QL.select_action(); // selecting an acting according to the initial policy
    QL.update_forcing(); // translating the action into physical parameters
    Fib.setforcing(QL.getp(), QL.getA()); //forcing the physical parameter
    if(learning == 1)
        QL.update_policy();
    
    while(it<nstep)
    {
        if((it % Nlearning)==0  && it != 0) {
            state = QL.compute_state(Fib.wind(U), Fib.orientation());
            QL.Qupdate(Fib.getcenter(0), previous_state);
            previous_state = state;
            if(learning == 1)
                QL.update_policy();
            QL.select_action();
            QL.update_forcing();
            Fib.setforcing(QL.getp(), QL.getA());
        }
        if((it % Nout)==0) {
            cout << "salam" << endl;
            cout<<setprecision(3);
            //cout << showpoint;
            cout<<"t = "<<t<<setw(10)<<"X = "<<Fib.getcenter(0)<<setw(10)<<"Vx = "<<Fib.getvelocity(0) << setprecision(1) << setw(10) << "Action: "<< QL.getaction()<<" State: "<< QL.getstate()<<endl;
            //Fib.save(outdir+"fiber"+i2s(it)+".nc",U);
        }
        if((it % Noutlearning)==0) {
            QL.save(it,t,outdir+"learn.bin");
        }
        Fib.evol(dt,U);
        t += dt;
        it++;
    }
    
    return 1;
}

