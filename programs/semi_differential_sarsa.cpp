//
//  learn2D_cellflow.cpp
//  i-eel
////  Created by Jeremie Bec on 24/02/2021.
//j

#include "src/Fiber2D.hpp"
#include "src/QLearning.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"
#include <armadillo>
#include "utils.h"

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
    const double Tmax = args.getreal("-T", "--time",1000000.0, "integration time");
    const int Nout = args.getint("-nout", "--step_out", 1000, "output period (in number of timesteps)");
    const double dt = args.getreal("-dt", "--timestep", 1e-3, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
    const double om = args.getreal("-om", "--frequency", 20, "Forcing frequency");
    const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
    const double u = args.getreal("-U", "--Velocity", 0.5, "Velocity amplitude");
    const double angle = args.getreal("-ang", "--angle", 0, "initial orientation: angle made with the horisontal in degrees");
    const int Nlearning = args.getint("-nl", "--step_learning", 10, "learning period (in number of timesteps)");
    const int Noutlearning = args.getint("-nlout", "--step_save_learning", 10, "learning output period (in number of timesteps)");
    const double gamma = args.getreal("-gamma", "--discountrate", 0.9995, "Discount rate");
    const double learnrate = args.getreal("-lr", "--learnrate", 0.005, "Learning rate");
    const double learnrate2 = args.getreal("-lr2", "--learnrate2", 0.001, "Learning rate");
    const double epsil = args.getreal("-eps", "--epsilon", 0.0, "Rate of random exploration");
    const double u0 = args.getreal("-slim", "--speed", 0.1, "Vitesse limite");
    //const double qinit = args.getreal("-q0", "--qinit", 0.25, "Initial Q entries");
    const int learning = args.getint("-lrn", "--learning", 1, "Input 1 for swimming with learning 0 otherwise");
    const int noflow = args.getint("-nfl", "--noflow", 0, "Input 1 for no flow 0 for cellular flow"); 
    const int incl_buckl = args.getint("-bckl", "--incl_buckl", 0, "Input 1 to include buckled states 0 otherwise");
    const int out = args.getint("-out", "--make_output", 0, "1 for saving trajectory files 0 otherwise");
    const int n_bootstrapping = args.getint("-nb", "--n_bootstrapping", 10000, "number of steps for bootstrapping");
    const int debug = args.getint("-debug", "--debug", 0, "1 for debug mode 0 otherwise");
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
    p.at(0) = -cos(angle * 3.14 / 180);
    p.at(1) = sin(angle * 3.14 / 180);
    cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
    Fiber2D Fib(Ns,L,zeta,E,beta,U,p);
    // set the forcing
    cout<<endl<<"------------------------------------------------"<<endl;
    cout<<"k="<<k<<"\tom="<<om<<"\talpha="<<alpha<<endl;
    Fib.setforcing(k, om);
    
    // define the learning
    MyMat Q;
    MyMat Pi;
    
    
    if(strcmp(Qinitdir.c_str(),"")==0){
        Q.load("input_data/Q.csv", arma::csv_ascii);
        cout<<"initial Q matrix read from csv file:" << endl;
        cout<<Q;
    }
    else{
        Q.load(Qinitdir, arma::csv_ascii);
    }

    if(strcmp(Qinitdir.c_str(),"")>0) {
        cout<<"Start Q from file "<<indir<<"learn.bin"<<endl;
        Q = readlastQ(indir+"learn.bin",12,8);
    }

    Pi.load("input_data/Pi.csv", arma::csv_ascii);
    cout<<"initial policy matrix read from csv file:" << endl;
    Q.zeros();
    Q = 2 * Q;
    //Pi.ones();
    cout<<Pi;
    

    MyCol Ampl;
    double a0 = zeta*alpha*om/(2.0*M_PI*(double)k/L);
    Ampl.set_size(4);
    Ampl(0) = 0;
    Ampl(1) = a0/3.0;
    Ampl(2) = 2.0*a0/3.0;
    Ampl(3) = a0;

    
    QLearning QL(Q, Pi, gamma,learnrate,u0,Ampl,epsil, true);
    char cname[512];
    string fname = outdir+"learn.bin";
    strcpy(cname, fname.c_str());
    FILE *fout = fopen(cname,"w");
    fclose(fout);
    

    Col<double> rewards; // a vector to store the last n rewards
    rewards.set_size(n_bootstrapping);
    rewards.zeros();
    Col<int> actions; // a vector to store the last n actions
    actions.set_size(n_bootstrapping);
    Col<int> states; // a vector to store the last n states
    states.set_size(n_bootstrapping);

    // time loop
    cout<<endl<<"------------------------------------------------"<<endl;
    int nstep = round(Tmax/dt);
    cout<<"Starting the time loop over "<<nstep<<" steps"<<endl;
    cout<<"output every "<<Nout<<" steps"<<endl;
    int it = 0;
    double t = 0;
    
    int state;
    int action;
    double reward;

    //int state = QL.compute_state(Fib.wind(U), Fib.orientation(), incl_buckl * Fib.calc_buckle());//initial state
    //QL.select_action();
    //int action = QL.get_action();
    //QL.update_forcing(); // translating the action into physical parameters
    //Fib.setforcing(QL.getp(), QL.getA()); //forcing the physical parameter
    Fib.save(outdir+"fiber"+i2s(0)+".ff",U);

    int learning_step; // a variable to store the step in multiples of Nlearn   
    double G; // to store the target value for TD(n) 
    double Rbar=0;
    double delta=0;
    int ind = 0;

    while(it<nstep)
    {
        if(((it % Nlearning))==0){
            learning_step = it / Nlearning;
            state = QL.compute_state(Fib.wind(U), Fib.orientation(), incl_buckl * Fib.calc_buckle()); // compute the new state
            QL.reward(Fib.getcenter(0));
            reward = QL.get_reward(); // reward of the previous step
            QL.select_action();
            action = QL.get_action();
            QL.update_forcing();
            Fib.setforcing(QL.getp(), QL.getA());
            
            //compute ind
            ind = learning_step % n_bootstrapping; // the index of the oldest stored values

            //learn
            if(learning_step >= n_bootstrapping){
                // update the target value
                G = 0;
                /*
                for(int i = 0; i < n_bootstrapping - 1; i++){
                    G = G + rewards[(ind + i + 1) % n_bootstrapping] - Rbar;
                }
                */

                G = sum(rewards) - n_bootstrapping * Rbar;
                // S_{t-n}, A_{t-n}
                int s = states.at(ind), a = actions.at(ind);
                //S_{t}, A_{t}
                int ind_prim = (ind + n_bootstrapping - 1) % n_bootstrapping;
                // old value of Q:
                double Q_old = QL.get_Q(s, a);
                //compute delta
                delta = G +  QL.get_Q(states[ind_prim], actions[ind_prim])  - Q_old;
                // update Rbar
                Rbar += learnrate2 * delta;
                // new estimate q_new = q_old + alpha * delta
                double Q_new = Q_old + learnrate * delta;
                if(it % Nout == 0) cout << "target: " << G + QL.get_Q(states[ind_prim], actions[ind_prim]) << " current estimate: " << Q_old << " G: " << G << " TD(n) error: " << delta << endl;
                QL.update_Q(s, a, Q_new);
                if(learning == 1) QL.update_policy();
                if (debug == 1) cout << " state: " << s << " action: " << a << " Q_old: " << Q_old << " G: " << G << " Q_new: " << Q_new << " Q new check: " << QL.get_Q(s, a) << endl;
            }

            //update stored values
            states.at(ind) = state;
            rewards.at(ind) = reward;
            actions.at(ind) = action;
            if(ind == 0)
            {
                QL.show_Q();
                QL.show_policy();
            }
        }
        

        
        if((it % Nout == 0)) {
            cout<<setprecision(3);
            //cout << showpoint;
            cout<<"t = "<<t<<setw(10)<<"X = "<<Fib.getcenter(0)<<setw(10)<<"Vx = "<<Fib.getvelocity(0) << setprecision(1) << setw(10) << "Action: "<< QL.getaction()<<" State: "<< QL.getstate()<< " mean perf: " << Rbar <<endl;
            cout << "size: " << rewards.n_rows << endl;
            if(out == 1) 
                {Fib.save(outdir+"fiber"+i2s(it)+".ff",U);}
        }
        if(((it) % Noutlearning)==0) {
            QL.save(it, t, outdir+"learn.bin");
            save(it, t, Fib.getcenter(0), state, action, reward, delta, Rbar, Q);
        }
        Fib.evol(dt,U);
        t += dt;
        it++;
    }
    return 1;
}