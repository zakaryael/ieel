//
//  learn2D_cellflow.cpp
//  i-eel
//
//  Created by Jeremie Bec on 24/02/2021.
//j


#include "src/Fiber2D.hpp"
#include "src/QLearning.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"

using namespace std;
using namespace arma;


MyCol get_features_vector(MyCol state, double x, double y){
    MyCol features; //4 state, 1 x, 1 y, 1 x % 2, 1 y % 2, total = 8 + 4 * 7 = 36
    int n = 8;
    features.set_size(n);
    MyCol features_with_interactions;
    features_with_interactions.set_size((std::floor(n * (n+1) / 2)));

    features.at(0) = state.at(0);
    features.at(1) = state.at(1);
    features.at(2) = state.at(2);
    features.at(3) = (float) state.at(3);
    features.at(4) = fmod(x, 2);
    features.at(5) = fmod(y , 2);
    features.at(6) = (float)fmod(x, 2) > 1;
    features.at(7) = (float)fmod(y, 2) > 1 ;

    for(int i=0; i < n; i++){
        features_with_interactions.at(i) = features.at(i);
        for(int j=i+1; j < n; j++){
            int ind = (j-i) * n + i - std::floor((j-i) * (j-i-1) / 2);
            features_with_interactions.at(ind) = features.at(i) * features.at(j);
        }
    }
    return features_with_interactions;
}

MyMat compute_policy(MyMat theta){
    MyMat P = exp(theta);
    for(int i=0; i < 6; i++){
        P.row(i) = P.row(i) / sum(P.row(i));
    }
    return P;
}

int select_action(MyMat policy, int state){
    //selects an action based on the policy
    time_t tt;
    std::default_random_engine rng;
    rng.seed((unsigned) time(&tt));
    arma::Row<double> tmp = policy.row(state);
    std::discrete_distribution<int> dist(tmp.begin(), tmp.end());
    return dist(rng);
}

int main(int argc, char* argv[]) {
    WriteProcessInfo(argc, argv);
    
    string purpose("Q-Learning oflearne\2D fiber in a cellular flow.");
    
    ArgList args(argc, argv, purpose);
    
    // Read the command-line options
    args.section("Program options");
    const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
    const double L = args.getreal("-L", "--length", 1.0, "fiber length");
    const double zeta = args.getreal("-z", "--zeta", 5e4, "friction coefficient");
    const double E = args.getreal("-E", "--EI", 1.0, "Young modulus");
    const double beta = args.getreal("-beta", "--penalisation", 400, "penalisation of extensibility");
    const double Tmax = args.getreal("-T", "--time",10000.0, "integration time");
    const int Nout = args.getint("-nout", "--step_out", 2000, "output period (in number of timesteps)");
    const double dt = args.getreal("-dt", "--timestep", 0.0005, "time step");
    const int Ns = args.getint("-ns", "--Ns", 200, "number of points in the fiber's discretization");
    const int k = args.getint("-k", "--wavenumber", 2, "Forcing wavenumber");
    const double om = args.getreal("-om", "--frequency", 2, "Forcing frequency");
    const double alpha = args.getreal("-alpha", "--alpha", 1, "Force amplitude");
    const double u = args.getreal("-U", "--Velocity", 0.05, "Velocity amplitude");
    const double u0 = args.getreal("-slim", "--speed", 0.01, "Vitesse limite");
    const int Nlearning = args.getint("-nl", "--step_learning", 2000, "learning period (in number of timesteps)");
    const int Noutlearning = args.getint("-nlout", "--step_save_learning", 2000, "learning output period (in number of timesteps)");
    const double beta0 = args.getreal("-beta0", "--learning_param2", 0.01, "Discount rate");
    const double alpha0 = args.getreal("-alpha0", "--learning_param1", 0.1, "Learning rate");
    //const int learning = args.getint("-lrn", "--learning", 1, "Input 1 for swimming with learning 0 otherwise");
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
    

    p.at(0) = -1;
    p.at(1) = 0;

    cout<<"initial orientation: p = ("<<p.at(0)<<","<<p.at(1)<<")"<<endl;
    Fiber2D Fib(Ns,L,zeta,E,beta,U,p);
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

    
    
    char cname[512];
    string fname = outdir+"naive.bin";
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
    int out = 0;


    //start of learning
    QLearning learner(u0, Ampl);

    MyCol v;
    v.set_size(6);
    MyMat theta;
    theta.set_size(6, 7);
    theta.zeros(); // initialize theta
    theta.load("Pi.csv", arma::csv_ascii);
    theta = 4 * theta;
    v.zeros(); // initialize v
    MyMat pi;
    pi.set_size(6, 7);
    cout << "hello!" << endl;
    pi = compute_policy(theta);
    

    
    double reward=0;

    int state = learner.compute_state(Fib.wind(U), Fib.orientation());
    int new_state;
    int action = select_action(pi, state);
    int new_action;
    double delta;
    double gamma = 0.9995;
    
    cout << action << endl;
    learner.set_action_to(action); // selecting an intial action  
    learner.update_forcing(); // translating the action into physical parameters
    Fib.setforcing(learner.getp(), learner.getA()); //forcing the physical parameter
    
    while(it<nstep) 
    {
        if((it % Nlearning)==0  && it != 0) {
            
            //observe s' and R
            new_state = learner.compute_state(Fib.wind(U), Fib.orientation());
            
            learner.reward(Fib.getcenter(0)); //weird?
            reward = learner.get_reward();
            

            //sample A'
            pi = compute_policy(theta);
            cout << pi;
            new_action = select_action(pi, new_state);
            
            //compute TD(0) error
            delta = reward + gamma * v.at(new_state) - v.at(state);
            
            //update theta parameters
            for(int b=0; b < 7; b++){
                theta(state, b) = theta(state, b) - alpha0 * delta * pi(state, b);
            }
            theta(state, action) = theta(state, action) + alpha0 * delta;
            
            //update v parameters
            v(state) = v(state) + beta0 * delta;
            
            //update action and state
            state = new_state;
            action = new_action;
            
            //take action A
            learner.set_action_to(action);
            learner.update_forcing();
            Fib.setforcing(learner.getp(), learner.getA());
         }
        if((it % Nout)==0) {
            cout<<setprecision(3);
            //cout << showpoint;
            cout<<"t = "<<t<<setw(4)<< " X = "<<Fib.getcenter(0)<<setw(4)<<" Vx = "<<Fib.getvelocity(0) << setprecision(1) << setw(4) << " Action: "<< learner.getaction()<<" State: "<< learner.getstate()<<endl;
            if(out == 1) Fib.save(outdir+"fiber"+i2s(it)+".ff",U);
        }
        if((it % Noutlearning)==0) {
            FILE *fout = fopen(cname,"a");
            double tmp = (double) state;
            fwrite(&tmp, sizeof(double), 1, fout);

            tmp = (double)action;
            fwrite(&tmp, sizeof(double), 1, fout);
            
            tmp = Fib.getcenter(0);
            fwrite(&tmp, sizeof(double), 1, fout);

            tmp = reward;
            fwrite(&tmp, sizeof(double), 1, fout);

            fwrite(&theta(0,0), sizeof(double), 42, fout);
            fwrite(&v(0,0), sizeof(double), 6, fout);
            cout << "__________________" << endl;
            fclose(fout);
        }
        Fib.evol(dt,U);
        t += dt;
        it++;
        
    }
    return 1;
}
