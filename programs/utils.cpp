
#include "utils.h"

void save(const std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar)
{
    char cname[512];
    strcpy(cname, file.c_str());
    FILE *fout = fopen(cname,"a");
    double tmp = (double)iteration;
    fwrite(&tmp, sizeof(double), 1, fout);
    double temp = time;
    fwrite(&temp, sizeof(double), 1, fout);
    temp = position;
    fwrite(&temp, sizeof(double), 1, fout);
    temp = (double) observation;
    fwrite(&temp, sizeof(double), 1, fout);
    temp = (double) action;
    fwrite(&temp, sizeof(double), 1, fout);
    temp = reward;
    fwrite(&temp, sizeof(double), 1, fout);
    temp = delta;
    fwrite(&temp, sizeof(double), 1, fout);
    temp = rbar;
    fwrite(&temp, sizeof(double), 1, fout);
    fclose(fout);
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
