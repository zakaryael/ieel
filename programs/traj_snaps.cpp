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
    const string Qinitdir = args.getpath("-Qdir", "--Qinitdir", "", "input directory to a csv file from which the initial Q is read");
    const string Pinitdir = args.getpath("-Pidir", "--Pinitdir", "", "input directory to a csv file from which the initial policy is read");  
    const double L = args.getreal("-L", "--length", 1.0, "fiber length");
    const double zeta = args.getreal("-z", "--zeta", 2.5e3, "friction coefficient");
    const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
    const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
    const double Tmax = args.getreal("-T", "--time",1.0e7, "integration time in multiples of nl");
    const double dt = args.getreal("-dt", "--timestep", 1e-3, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
    const double om = args.getreal("-om", "--frequency", 20, "Forcing frequency");
    const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
    const double u = args.getreal("-U", "--Velocity", 0.5, "Velocity amplitude");
    const int noflow = args.getint("-nfl", "--noflow", 0, "Input 1 for no flow 0 for cellular flow"); 
    const int incl_buckl = args.getint("-bckl", "--incl_buckl", 0, "Input 1 to include buckled states 0 otherwise");
    const double u0 = args.getreal("-slim", "--speed", 0.1, "Vitesse limite");
    const int Nlearning = args.getint("-nl", "--step_learning", 100, "learning period (in number of timesteps)");

    args.check();
    mkdir(outdir);
    args.save(outdir);
    
    // fluid flow set to cellular with amplitude u
    Flow2D U(Cellular);
    U.initcellular(u);
    // To set the fluid flow set to 0:
    if (noflow == 1){Flow2D U(Null);}
    


    Fiber2D Fib(Ns, L, zeta, E, beta);
    // set the forcing
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
    Fib.setforcing(k, om);
    
    // define the learning

    MyCol Ampl;
    double a0 = zeta*alpha*om/(2.0*M_PI*(double)k/L);
    Ampl.set_size(4);
    Ampl(0) = 0;
    Ampl(1) = a0/3.0;
    Ampl(2) = 2.0*a0/3.0;
    Ampl(3) = a0;

    
    QLearning QL(u0,Ampl, 7);
    char cname[512];
    string fname = outdir+"naive_snaps.bin";
    strcpy(cname, fname.c_str());
    FILE *fout = fopen(cname,"w");
    fclose(fout);
    
    // time loop
    cout<<endl<<"------------------------------------------------"<<endl;
    //int nstep = round(1/dt) / ;
    cout<<"output every 1 sec"<<endl;
    int it = 0;

    for (int n=0; n < Tmax; n++) 
    {  
        for(int a=0; a < 7; a++){
        FILE *fout = fopen(cname,"a");
        // Appends 6 doubles to the file: time, position, current_state, action, next_state, reward.
        double tmp = (double)n;
        fwrite(&tmp, sizeof(double), 1, fout);
        QL.set_action_to(a);

        tmp = (double)a;
        fwrite(&tmp, sizeof(double), 1, fout);
        QL.update_forcing();
        Fib.read(U, indir+"fiber"+i2s(n * Nlearning)+".ff",QL.getp(),QL.getA(), n);
        
        int current_state = QL.compute_state(Fib.wind(U), Fib.orientation(), incl_buckl * Fib.calc_buckle());
        cout<<"t = "<<n<<setw(10)<<"X = " << setprecision(4)<<Fib.getcenter(0)<<setw(10)<<"Vx = "<<Fib.getvelocity(0) << setprecision(1) << setw(10) << " State: "<< current_state<<endl;
        cout<<"Action: "<< QL.getaction()<<endl;
        double xold = Fib.getcenter(0);
        fwrite(&xold, sizeof(double), 1, fout);

        tmp = (double)current_state;
        fwrite(&tmp, sizeof(double), 1, fout);

        it = 0;
        while(it <= Nlearning){
            Fib.evol(dt,U);
            it++;
        }
        int next_state = QL.compute_state(Fib.wind(U), Fib.orientation(), incl_buckl * Fib.calc_buckle()); 
        cout<<"t = "<<n+1<<setw(10)<<"X = "<< setprecision(4)<<Fib.getcenter(0)<<setw(10)<<"Vx = "<<Fib.getvelocity(0) << setprecision(1) << setw(10) << " State: "<< next_state<<endl; 
        
        tmp = (double)next_state;
        fwrite(&tmp, sizeof(double), 1, fout);

        tmp = Fib.getcenter(0);
        fwrite(&tmp, sizeof(double), 1, fout);
        tmp = tmp - xold;
        cout << "reward: " << tmp << endl;
        fwrite(&tmp, sizeof(double), 1, fout);
        fclose(fout);
        cout << "__________________" << endl;
        }
        cout << "------------------" << endl << endl;
    }
    return 1;
}
