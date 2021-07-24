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


QLearning::QLearning(MyMat Q, double gamma, double learnrate, double u0, MyCol Ampl, double epsil){
    if(Q.n_cols != 2*Ampl.n_rows)
        ErrorMsg("initQlearning: Q and Ampl are incompatible");
    Q_ = Q;
    nstate_ = Q.n_rows;
    naction_ = Q.n_cols;
    gamma_ = gamma;
    learnrate_ = learnrate;
    u0_ = u0;
    Ampl_ = Ampl;
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

double QLearning::state_update(double u, vector<double> P, int buckl) {
    // Compute the new state
    state_ = discr_wind(u)+3*(buckl+2*discr_orientation(P));
    // Maximize Q and find the action
    double qmax = -10;
    for(int ia=0; ia<naction_; ++ia) {
        if(Q_(state_,ia) > qmax) {
            qmax = Q_(state_,ia);
            action_ = ia;
        }
    }
    
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    std::uniform_int_distribution<> int_unif(0, 2*Ampl_.n_rows-1);
    float v = uniform(generator_);
    int action_new;
    // std::cout << "\n the random picked number is: "<< v<< endl;
    if(v < epsilon_) {// if v < epsilon pick the action at random
        action_new = int_unif(generator_);
        while(action_new==action_)
            action_new = int_unif(generator_);
        action_ = action_new;
    }
    // Update the forcing parameters depending on the action
    if(action_<(int)Ampl_.n_rows) {
        p_.at(0)=0; p_.at(1)=1;
    }
    else {
        p_.at(0)=1; p_.at(1)=0;
    }
    A_ = Ampl_(action_%Ampl_.n_rows);
    return qmax;
}


void QLearning::Qupdate(double xnew, double u, vector<double> P, int buckl) {
    int stateold = state_;
    int actold = action_;
    double qmax = state_update(u,P,buckl);
    reward(xnew);
    Q_(stateold,actold) = (1.0-learnrate_)*Q_(stateold,actold)+learnrate_*(rew_+gamma_*qmax);
}

void QLearning::save(int istep, const std::string& filebase) {
    char cname[512];
    strcpy(cname, filebase.c_str());
    FILE *fout = fopen(cname,"a");
    // Appends 5 + nstate_*naction_ doubles to the file
    double tmp = (double)istep;
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

