//
//  QLearning.hpp
//  i-eel
//
//  Created by Jeremie Bec on 23/07/2021.
//

#ifndef QLearning_hpp
#define QLearning_hpp
#include "basics/Utilities.h"

using namespace std;
using namespace arma;

typedef Mat<double> MyMat;
typedef Col<double> MyCol;

MyMat readlastQ(const std::string& filebase, int ns, int na);

class QLearning {
public:
    // Class initialisation
    QLearning();
    QLearning(double u0, MyCol Ampl);
    QLearning(MyMat Q, MyMat, double gamma, double learnrate, double u0, MyCol Ampl, double epsil=0.0, bool merge_zeros=false);
    // 
    void set_seed(void);
    int discr_wind(double wind);
    int discr_orientation(vector<double> P);
    void save(const std::string& filebase);
    
    //
    int compute_state(double u, vector<double> P, int buckl);
    int compute_state(double u, vector<double> P);
    void select_action(void);
    int get_action(void);
    void update_forcing(void);
    inline void reward(double xnew){
        rew_ = xnew-xold_;
        xold_ = xnew;
    };

    double get_reward(){
        return rew_;
    };
    void show_Q(void){
        std::cout << "Q: " << endl;
        std::cout << Q_;
    };
    void show_policy(void){
        std::cout << "pi: " << endl;
        std::cout << policy_;
    };
    double get_Q(int, int); // get the value of the q matrix for a state-action pair
    // update the Q matrix
    void update_Q(int, int, double); // updates the Q value corresponding to a state-action pair
    void Qupdate(double xnew, int);
    void update_policy(void);
    // define the reward



    inline double printreward(){ return rew_; };
    inline vector<double> getp() { return p_; };
    inline double getA() { return A_; };
    inline double getaction() { return action_; };
    inline double getstate() { return state_; };

    void set_action_to(int a);
    void save(int istep, double time, const std::string& filebase);
    
private:
    // Learning parameters
    MyMat Q_; // Q-table
    MyMat policy_;
    int nstate_; // number of states
    int naction_; // number of actions
    int state_;
    int action_;
    double rew_;
    double gamma_; // discount rate
    double learnrate_; // learning rate
    double u0_; // discretization of the wind
    float epsilon_ = 0.0;  // exploration rate
    std::default_random_engine generator_;
    MyCol Ampl_; // table containing the discrete values of the forcing amplitude
    vector<double> p_; // direction of the helix
    double A_ = 0.0;   // internal amplitude
    double xold_ = 0; // needed to compute the reward
    bool merge_zeros_ = true;
};

#endif /* QLearning_hpp */