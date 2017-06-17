#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
public:
    vector<double> x_predicted;
    vector<double> y_predicted;

    MPC();

    virtual ~MPC();

    // Retrieve Lf.
    double getLf();

    // Initialize x_predicted and y_predicted.
    void Init();

    // Solve the model given an initial state and polynomial coefficients.
    // Return the first actuations.
    vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
