#include "src/Fiber2D.hpp"
#include "src/QLearning.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"
#include <sys/stat.h>


using namespace std;

int main(int argc, char* argv[]) {
    WriteProcessInfo(argc, argv);
    
    string purpose("Q-Learning of a 2D fiber in a cellular flow.");
    
    ArgList args(argc, argv, purpose);
    
    // Read the command-line options
    args.section("Program options");
    const string wdir = args.getpath("-dir", "--directory", "output_data/wdir/", "work directory");
    const int iteration = args.getint("-it", "--iteration", 0, "iteration tracker");
    const double L = args.getreal("-L", "--length", 1.0, "fiber length");
    const double zeta = args.getreal("-z", "--zeta", 2.5e3, "friction coefficient");
    const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
    const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
    //const double Tmax = args.getreal("-T", "--time",1.0e7, "integration time in multiples of nl");
    const double dt = args.getreal("-dt", "--timestep", 1e-3, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
    const double om = args.getreal("-om", "--frequency", 20, "Forcing frequency");
    const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
    const double u = args.getreal("-U", "--Velocity", 0.5, "Velocity amplitude");
    const int noflow = args.getint("-nfl", "--noflow", 0, "Input 1 for no flow 0 for cellular flow"); 
    //const int incl_buckl = args.getint("-bckl", "--incl_buckl", 0, "Input 1 to include buckled states 0 otherwise");
    const double u0 = args.getreal("-slim", "--speed", 0.1, "Vitesse limite");
    const int Nswim = args.getint("-nsw", "--step_swimming", 100, "swimming period (in number of timesteps)");
	const int action = args.getint("-a", "action", 6, "the action to implement by the swimmer");
    const int nout = args.getint("-nout", "--step_out", 100, "output period to command line(in number of timesteps)");

    args.check();

    /* checks if a dir exists
    struct stat sb;
    if (stat(wdir.c_str(), &sb) == 0) //&& S_ISDIR(sb.st_mode)
    {
        cout << ":)" << endl;
    }
    */

    // fluid flow set to cellular with amplitude u
    Flow2D U(Cellular);
    U.initcellular(u);
    // To set the fluid flow set to 0:
    if (noflow == 1){Flow2D U(Null);}

    // amplitude
    MyCol Ampl;
    double a0 = zeta*alpha*om/(2.0*M_PI*(double)k/L);
    Ampl.set_size(4);
    Ampl(0) = 0;
    Ampl(1) = a0/3.0;
    Ampl(2) = 2.0*a0/3.0;
    Ampl(3) = a0;

    //deal with the initial iteration
    Fiber2D Fib;
    mkdir(wdir);
    if(iteration == 0){
        args.save(wdir);
        vector<double> p(2);
        p.at(0) = -1;
        p.at(1) = 0;
        Fib = Fiber2D(Ns,L,zeta,E,beta,U,p);
        Fib.save(wdir+"fiber"+i2s(iteration)+".ff",U);
        cout<<endl<<"------------------------------------------------"<<endl;
        cout<< "files are stored in: " << wdir.c_str() << endl;
        cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
        cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
        

    }
    else{
        Fib = Fiber2D(Ns, L, zeta, E, beta);
        Fib.read(U, wdir+"fiber"+i2s(iteration)+".ff");
    }
    
    // set the forcing
    Fib.setforcing(k, om);

    //create the swimmer
    QLearning Swimmer(u0,Ampl);
    
    //set and implement the action
    Swimmer.set_action_to(action);
    Swimmer.update_forcing();
    Fib.setforcing(Swimmer.getp(), Swimmer.getA());
    

    for(int it = 0; it < Nswim; it++){
        Fib.evol(dt,U);
    }
    Fib.save(wdir+"fiber"+i2s(iteration+1)+".ff",U);

    if(iteration % nout == 0){cout << "iteration " << iteration << endl;}
    return 1;
}
