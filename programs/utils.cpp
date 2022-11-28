
#include "utils.h"

void save(const std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar)
{
    char cname[512];
    strcpy(cname, file.c_str());
    FILE *fout = fopen(cname,"a");
    // Appends ... to the file
    double tmp = (double)iteration;
    fwrite(&tmp, sizeof(double), 1, fout);
    fwrite(&time, sizeof(double), 1, fout);
    fwrite(&position, sizeof(double), 1, fout);
    fwrite(&observation, sizeof(double), 1, fout);
    fwrite(&action, sizeof(double), 1, fout);
    fwrite(&reward, sizeof(double), 1, fout);
    fwrite(&delta, sizeof(double), 1, fout);
    fwrite(&rbar, sizeof(double), 1, fout);
    fclose(fout);

}

void save(const std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar, MyMat Q){
    int n_act = Q.n_cols;
    int n_obs = Q.n_rows;
    save(file, iteration, time, position, observation, action, reward, delta, rbar);
    char cname[512];
    strcpy(cname, file.c_str());
    FILE *fout = fopen(cname,"a");
    
    fwrite(&Q(0,0), sizeof(double), n_act * n_obs, fout);
    fclose(fout);

}