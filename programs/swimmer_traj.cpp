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
   // const int iteration = args.getint("-it", "--iteration", 0, "iteration tracker");
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
    const double u0 = args.getreal("-slim", "--speed", 0.1, "Vitesse limite");
    const int Nswim = args.getint("-nsw", "--step_swimming", 10, "swimming period (in number of timesteps)");
    const int nout = args.getint("-nout", "--step_out", 1, "output period to command line(in number of timesteps)");
    const unsigned int Tmax = args.getint("-T", "--Tmax", 1000000000, "final time");


    args.check();
    double x_old;


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

    //creating a binary file where the data will go
    char cname[512];
    string fname = wdir+"snaps.bin";
    strcpy(cname, fname.c_str());
    FILE *fout = fopen(cname,"w");
    fclose(fout);
    Fiber2D Fib(Ns, L, zeta, E, beta);//init the fiber
    int state;
    // set the forcing
    Fib.setforcing(k, om);
    //create the swimmer
    QLearning Swimmer(u0,Ampl);
    for(int iteration=0; iteration < Tmax * Nswim; iteration++)
        {
        for(int action=0; action < 7; action++){
            //set and implement the action
            Swimmer.set_action_to(action);
            Swimmer.update_forcing();//updates p_ and A_
            Fib.setforcing(Swimmer.getp(), Swimmer.getA());
            Fib.read(U, wdir+"fiber"+i2s(iteration*Nswim)+".ff", iteration*Nswim*dt, dt);
            state = Swimmer.compute_state(Fib.wind(U), Fib.orientation(),0);
            x_old = Fib.getcenter(0);

            if(iteration % nout == 0){cout << "position of the center at iteration " << iteration << " (t= " << (iteration) * Nswim * dt << "s) is: " << Fib.getcenter(0) << endl;}
            
            for(int it = 0; it < Nswim; it++){
                Fib.evol(dt,U);
            }
            
            // saving the data
                //one file per iteration
            //Fib.save(wdir+"fiber"+i2s(iteration+1)+"action"+i2s(action)+".ff",U);
                //same file for the whole run:
                
                if((state <= 3 && action == 3) || (state >=4 && action == 6))
                {
                FILE *fout = fopen(cname,"a");
                double tmp = (double)iteration;
                fwrite(&tmp, sizeof(double), 1, fout); // iteration
                tmp = (double)iteration * Nswim * dt;
                fwrite(&tmp, sizeof(double), 1, fout); //time
                tmp = (double)action;
                fwrite(&tmp, sizeof(double), 1, fout); //the action
                tmp = Fib.getcenter(0);
                fwrite(&tmp, sizeof(double), 1, fout); //the position
                fwrite(&x_old, sizeof(double), 1, fout);//the old position
                fclose(fout);
                }
        }
    }
    return 1;
}
