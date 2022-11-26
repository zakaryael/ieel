//
//  QLearning.cpp
//  i-eel
//
//  Created by Jeremie Bec on 23/07/2021.
//

#include "QLearning.hpp"

MyMat readlastQ(const std::string& filebase, int ns, int na) {
    MyMat Q(ns,na);
    int it=0, sz=ns*na, rdsz;
    double tmp[5];
    char cname[512];
    strcpy(cname, filebase.c_str());
    FILE *fin = fopen(cname,"r");
    fread(&tmp[0],sizeof(double),5,fin);
    while(!feof(fin)) {
        rdsz = fread(&Q[0],sizeof(double),sz,fin);
        fread(&tmp[0],sizeof(double),5,fin);
        it++;
    }
    if(rdsz!=sz)
        ErrorMsg("it = "+i2s(it)+": Could not read from file "+filebase);
    cout<<"Detected "<<it<<" iterations"<<endl;
    fclose(fin);
    return Q;
}

QLearning::QLearning(double u0, MyCol Ampl){
    u0_ = u0;
    Ampl_ = Ampl;
    p_.resize(2);
}

QLearning::QLearning(MyMat Q, MyMat Pi, double gamma, double learnrate, double u0, MyCol Ampl, double epsil, bool merge_zeros){
    merge_zeros_ = merge_zeros;
    if(merge_zeros_) {
        if(Q.n_cols != 2*Ampl.n_rows-1){
            std::cout << "initQlearning: Q  ( " << Q.n_cols << " ) and Ampl ( " << 2 * Ampl.n_rows - 1<< " ) cols are incompatible" << endl;
            ErrorMsg("initQlearning: Q and Ampl rows are incompatible (merge)");
        }
   }
    else {
        if(Q.n_cols != 2*Ampl.n_rows){
            cout << "initQlearning: Q  ( " << Q.n_cols << " ) and Ampl ( " << Ampl.n_rows << " ) cols are incompatible" << endl;
            ErrorMsg("initQlearning: Q and Ampl rows are incompatible");
        }

    }
    Q_ = Q;
    nstate_ = Q.n_rows;
    naction_ = Q.n_cols;
    gamma_ = gamma;
    learnrate_ = learnrate;
    u0_ = u0; 
    Ampl_ = Ampl; // ??
    policy_ = Pi;
    epsilon_ = epsil;
    set_seed();
    p_.resize(2);
    
}

void QLearning::set_seed(void){
    time_t tt;
    generator_.seed((unsigned) time(&tt));
}

int QLearning::discr_wind(double wind) {
    int iwind = 0;
    if(wind>-u0_) {
        if(wind<u0_)
            iwind = 1;
        else
            iwind = 2;
    }
    return iwind;
}

int QLearning::discr_orientation(vector<double> P){
    int iorient = 0;
    if(P.at(0)>0)
        iorient = 1;
    return iorient;
}

int QLearning::compute_state(double u, vector<double> P, int buckl){
    // Computes and updates the current state
    int state = discr_wind(u)+3*(2 * buckl + discr_orientation(P));
    state_ = state;
    return state;
}

int QLearning::compute_state(double u, vector<double> P){
    // Computes and updates the current state
    int state = 2*discr_orientation(P);
    if(discr_wind(u)>0)
        state++;
    state_ = state;
    return state;
}

void QLearning::select_action(void){
    //selects an action based on the currently followed policy
    arma::Row<double> tmp = policy_.row(state_);
    std::discrete_distribution<int> dist(tmp.begin(), tmp.end());
    action_ = dist(generator_);
}
int QLearning::get_action(void){
    return action_;
}

void QLearning::set_action_to(int a){
    action_ = a;
}
void QLearning::update_policy(void){
    //updates the followed policy
    //epsilon-greedy policy 
    policy_.row(state_).fill(epsilon_ / (double)(naction_ - 1));
    policy_(state_, Q_.row(state_).index_max()) = 1 - epsilon_;
}

void QLearning::update_Q(int state, int action, double value){
    Q_.at(state, action) = value;
}

double QLearning::get_Q(int state, int action){
    return Q_.at(state, action);
}

void QLearning::Qupdate(double xnew, int previous_state) {
    //updates the value of Q
    reward(xnew); // :(
    Q_(previous_state, action_) = (1.0 - learnrate_) * Q_(previous_state, action_) + learnrate_ * (rew_ + gamma_ * Q_.row(state_).max());
}

void QLearning::update_forcing(void){
    // Update the forcing parameters depending on the action
    if(merge_zeros_) {
        naction_ = 2 * Ampl_.size() - 1;
        int ii = action_ -(naction_-1)/2;
        
        if(ii>=0) {
            p_.at(0)=1; p_.at(1)=0;
            A_ = Ampl_(ii);
        }
        else {
            p_.at(0)=0; p_.at(1)=1;
            A_ = Ampl_(-ii);
        }
    }
    else {
        if(action_<(int)Ampl_.n_rows) {
            p_.at(0)=0; p_.at(1)=1;
        }
        else {
            p_.at(0)=1; p_.at(1)=0;
        }
        A_ = Ampl_(action_%Ampl_.n_rows);
    }
}

void QLearning::save(int istep, double time, const std::string& filebase) {
    char cname[512];
    strcpy(cname, filebase.c_str());
    FILE *fout = fopen(cname,"a");
    // Appends 5 + nstate_*naction_ doubles to the file
    double tmp = (double)istep;
    fwrite(&time, sizeof(double), 1, fout);
    fwrite(&tmp, sizeof(double), 1, fout);
    fwrite(&xold_, sizeof(double), 1, fout);
    fwrite(&rew_, sizeof(double), 1, fout);
    tmp = (double)state_;
    fwrite(&tmp, sizeof(double), 1, fout);
    tmp = (double)action_;
    fwrite(&tmp, sizeof(double), 1, fout);
    fwrite(&Q_(0,0), sizeof(double), nstate_*naction_, fout);
    fclose(fout);
}

