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
    const string outdir = args.getpath("-o", "--outdir", "output_data/untiteled/", "output directory");
    const string indir = args.getpath("-in", "--indir", "", "input directory to a binary file from which initial Q is read");
    const string PiFile = args.getpath("-Pi", "--PiFile", "", "csv file from which the initial policy is read");
    const double L = args.getreal("-L", "--length", 1.0, "fiber length");
    const double zeta = args.getreal("-z", "--zeta", 2.5e3, "friction coefficient");
    const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
    const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
    const double dt = args.getreal("-dt", "--timestep", 1e-3, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
    const double om = args.getreal("-om", "--frequency", 20, "Forcing frequency");
    const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
    const double u = args.getreal("-U", "--Velocity", 0.5, "Velocity amplitude");
    const int noflow = args.getint("-nfl", "--noflow", 0, "Input 1 for no flow 0 for cellular flow");
    const int incl_buckl = args.getint("-bckl", "--incl_buckl", 0, "Input 1 to include buckled states 0 otherwise");
    const double u0 = args.getreal("-slim", "--speed", 0.1, "Vitesse limite");
    const int Nlearn = args.getint("-nl", "--step_learning", 100, "learning period (in number of timesteps)");
    const int Ntot = args.getint("-Ntot", "--num_step",10000, "integration time in units of Nlearn*dt");
    
    args.check();
    mkdir(outdir);
    args.save(outdir);
    
    // fluid flow set to cellular with amplitude u
    Flow2D U(Cellular);
    U.initcellular(u);
    // To set the fluid flow set to 0:
    if (noflow == 1){Flow2D U(Null);}
    
    
    
    // define the learning
    vector<int> PiVec = {3,3,3,3,6,6};
    vector<double> pinit = {-1,0};
    Fiber2D Fib(Ns, L, zeta, E, beta, U, pinit);
    vector <Fiber2D> FibTab(7);
    // set the forcing
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
    Fib.setforcing(k, om);
    
    
    MyCol Ampl;
    double a0 = zeta*alpha*om/(2.0*M_PI*(double)k/L);
    Ampl.set_size(4);
    Ampl(0) = 0;
    Ampl(1) = a0/3.0;
    Ampl(2) = 2.0*a0/3.0;
    Ampl(3) = a0;

    
    QLearning QL(u0,Ampl);
    char cname[512];
    string fname = outdir+"naive_snaps.bin";
    strcpy(cname, fname.c_str());
    FILE *fout = fopen(cname,"w");
    fclose(fout);
    
    // time loop
    cout<<endl<<"------------------------------------------------"<<endl;
    for (int n=0; n < Ntot; n++) 
    {
        double xold = Fib.getcenter(0);
        int oldstate = QL.compute_state(Fib.wind(U), Fib.orientation(), incl_buckl * Fib.calc_buckle());
        
        for(int a=0; a<7; a++) {
            FibTab.at(a) = Fib;
            QL.set_action_to(a);
            QL.update_forcing();
            FibTab.at(a).setforcing(QL.getp(), QL.getA());
            //Fib.setforcing(QL.getp(), QL.getA());
            for(int it=0; it<Nlearn; ++it) {
                //cout << "entered the loop!" << endl;
                //Fib.evol(dt, U);
                FibTab.at(a).evol(dt,U);
                //cout << "passed!" << endl;
            }
            // Outputs
            FILE *fout = fopen(cname,"a");
            // 0: step number
            double tmp = (double)n;
            fwrite(&tmp, sizeof(double), 1, fout);
            // 1: time
            tmp = (double)(n*Nlearn)*dt;
            fwrite(&tmp, sizeof(double), 1, fout);
            // 2: action
            tmp = (double)a;
            fwrite(&tmp, sizeof(double), 1, fout);
            // 3: former position of the center of mass
            fwrite(&xold, sizeof(double), 1, fout);
            // 4: former state
            tmp = (double)oldstate;
            fwrite(&tmp, sizeof(double), 1, fout);
            // 5: new state
            tmp = (double)QL.compute_state(FibTab.at(a).wind(U), FibTab.at(a).orientation(), incl_buckl * FibTab.at(a).calc_buckle());
            fwrite(&tmp, sizeof(double), 1, fout);
            // 6: new position of the center of mass
            tmp = FibTab.at(a).getcenter(0);
            fwrite(&tmp, sizeof(double), 1, fout);
            // 7: reward
            tmp = tmp - xold;
            fwrite(&tmp, sizeof(double), 1, fout);
            fclose(fout);
        }

        Fib = FibTab.at(PiVec.at(oldstate));
    }
cout << "__________________" << endl;
return 1;
}
