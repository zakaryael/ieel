#include "src/Fiber2D.hpp"
#include "src/QLearning.hpp"
#include "basics/RunsIO.h"
#include "basics/Arglist.h"
#include <armadillo>

void save(const std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar);
// void save(const std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar, MyMat Q);
// void save(std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar, MyMat Q, MyMat Pi);
// void save(std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar, MyCol v, MyMat Pi);
// void save(std::string& file, int iteration, double time, double position, double observation, double action, double reward, double delta, double rbar, MyCol v, MyMat Q, MyMat Pi);

MyMat compute_policy(MyMat theta);
int select_action(MyMat policy, int state);
