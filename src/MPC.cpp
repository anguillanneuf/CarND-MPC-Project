#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the time step length and duration
size_t N;
double dt;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

double ref_v;

// (same as in mpc_to_line quiz)
// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lives easier.
size_t x_start;
size_t y_start;
size_t psi_start;
size_t v_start;
size_t cte_start;
size_t epsi_start;
size_t delta_start;
size_t a_start;

// Adjust weight for each term of the cost function to get desired results.
double cte_cost_weight;
double epsi_cost_weight;
double speed_cost_weight;
double throttle_cost_weight;
double steering_cost_weight;
double delta_throttle_cost_weight;
double delta_steering_cost_weight;

class FG_eval {
public:
    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;

    FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector &fg, const ADvector &vars) {
        // TODO: implement MPC
        // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
        // NOTE: You'll probably go back and forth between this function and
        // the Solver function below.

        // The cost is stored is the first element of `fg`.
        // Any additions to the cost should be added to `fg[0]`.
        fg[0] = 0;

        // Reference State Cost
        // TODO: Define the cost related the reference state and anything else beneficial.

        // Consideration 1: desired position
        for (int t = cte_start; t < N + cte_start; t++) {
            fg[0] += cte_cost_weight * CppAD::pow(vars[t], 2);
        }
        // Consideration 2: desired heading
        for (int t = epsi_start; t < N + epsi_start; t++) {
            fg[0] += epsi_cost_weight * CppAD::pow(vars[t], 2);
        }

        /* adjust ref_v according to curve and steering angle
         * http://www.intmath.com/applications-differentiation/8-radius-curvature.php
         * */

        AD<double> dfdx = coeffs[1] + 2 * coeffs[2] * vars[x_start] + 3 * coeffs[3] * pow(vars[x_start], 2);
        AD<double> dfdx2 = 2 * coeffs[2] + 2 * 3 * coeffs[3] * vars[x_start];
        AD<double> R_curve = pow(1 + dfdx * dfdx, 1.5) / abs(dfdx2);

        if (vars[delta_start] > 4.0 * M_PI /180.0 || vars[delta_start] < -4.0 * M_PI /180.0 || R_curve < 50.0){
            ref_v = 55.0;
        } else {
            ref_v = 105.0;
        }

        // Consideration 3: dealing with stopping
        for (int t = v_start; t < N + v_start; t++) {
            fg[0] += speed_cost_weight * CppAD::pow(vars[t] - ref_v, 2);
        }
        // Considerations 4: add steering control input to avoid jerking the steering wheel
        for (int t = delta_start; t < N - 1 + delta_start; t++) {
            fg[0] += steering_cost_weight * CppAD::pow(vars[t], 2);
        }
        // Considerations 5: add throttle control input to avoid sudden change in acceleration
        for (int t = a_start; t < N - 1 + a_start; t++) {
            fg[0] += throttle_cost_weight * CppAD::pow(vars[t], 2);
        }
        // Consideration 6: adds diff between the next steering actuator state & the current one
        for (int t = delta_start; t < N - 2 + delta_start; t++) {
            fg[0] += delta_steering_cost_weight * CppAD::pow(vars[t+1] - vars[t], 2);
        }
        // Consideration 7: adds diff between the next acceleration actuator state & the current one
        for (int t = a_start; t < N - 2 + a_start; t++) {
            fg[0] += delta_throttle_cost_weight * CppAD::pow(vars[t+1] - vars[t], 2);
        }

        //
        // Setup Constraints
        //

        // Initial constraints
        //
        // We add 1 to each of the starting indices due to cost being located at
        // index 0 of `fg`.
        // This bumps up the position of all the other values.

        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + psi_start] = vars[psi_start];
        fg[1 + v_start] = vars[v_start];
        fg[1 + cte_start] = vars[cte_start];
        fg[1 + epsi_start] = vars[epsi_start];

        // The rest of the constraints
        for (int t = 1; t < N; t++) {
            AD<double> x1 = vars[x_start + t];
            AD<double> x0 = vars[x_start + t - 1];
            AD<double> psi0 = vars[psi_start + t - 1];
            AD<double> v0 = vars[v_start + t - 1];

            // `x`
            fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);

            // `y`
            AD<double> y1 = vars[y_start + t];
            AD<double> y0 = vars[y_start + t - 1];
            fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);

            // `psi`
            AD<double> psi1 = vars[psi_start + t];
            AD<double> delta0 = vars[delta_start + t - 1];
            fg[1 + psi_start + t] = psi1 - (psi0 + v0 / Lf * delta0 * dt);

            // `v`
            AD<double> v1 = vars[v_start + t];
            AD<double> a0 = vars[a_start + t - 1];
            fg[1 + v_start + t] = v1 - (v0 + a0 * dt);

            // `cte`
            AD<double> cte1 = vars[cte_start + t];
            AD<double> epsi0 = vars[epsi_start + t - 1];
            AD<double> cte0 = vars[cte_start + t - 1];
            AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
            fg[1 + cte_start + t] = cte1 - ((f0 - y0) + v0 * CppAD::sin(epsi0) * dt);

            // `epsi`
            AD<double> epsi1 = vars[epsi_start + t];
            AD<double> psides0 = CppAD::atan(coeffs[1] + (2 * coeffs[2] * x0) + (3 * coeffs[3] * CppAD::pow(x0, 2)));
            fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 / Lf * delta0 * dt);

        }
    }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}

MPC::~MPC() {}

double MPC::getLf(){
    return Lf;
};

void MPC::Init(double w0, double w1, double w2, double w3, double w4, double w5, double w6,
               int N_input, double dt_input, double v_input){

    // read from CL
    cte_cost_weight = w0;
    epsi_cost_weight = w1;
    speed_cost_weight = w2;
    throttle_cost_weight = w3;
    steering_cost_weight = w4;
    delta_throttle_cost_weight = w5;
    delta_steering_cost_weight = w6;
    N = N_input;
    dt = dt_input;
    ref_v = v_input;

    // update indices
    x_start = 0;
    y_start = x_start + N;
    psi_start = y_start + N;
    v_start = psi_start + N;
    cte_start = v_start + N;
    epsi_start = cte_start + N;
    delta_start = epsi_start + N;
    a_start = delta_start + N - 1;

    for(int i = 0; i < N - 1; i++){
        x_predicted.push_back(0.0);
        y_predicted.push_back(0.0);
    }
}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
    bool ok = true;
    size_t i;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    // TODO: Set the number of model variables (includes both states and inputs).
    // For example: If the state is a 4 element vector, the actuators is a 2
    // element vector and there are 10 timesteps. The number of variables is:
    //
    // 4 * 10 + 2 * 9
    size_t n_vars = N * 6 + (N - 1) * 2;

    // Number of constraints
    size_t n_constraints = N * 6;

    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];

    // Initial value of the independent variables.
    // SHOULD BE 0 besides initial state.
    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++) {
        vars[i] = 0;
    }

    // Set the initial variable values
    vars[x_start] = x;
    vars[y_start] = y;
    vars[psi_start] = psi;
    vars[v_start] = v;
    vars[cte_start] = cte;
    vars[epsi_start] = epsi;

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    for (int i = 0; i < delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // The upper and lower limits of delta are set to -25 and 25
    // degrees (values in radians).
    // NOTE: Feel free to change this to something else.
    for (int i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = -25.0 * M_PI /180.0; // vars in radians
        vars_upperbound[i] = 25.0 * M_PI /180.0;
    }

    // Acceleration/decceleration upper and lower limits.
    // NOTE: Feel free to change this to something else.
    for (int i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    // Set lower and upper limits for variables.

    // Lower and upper limits for the constraints
    // Should be 0 besides initial state.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (int i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }

    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[epsi_start] = epsi;

    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[epsi_start] = epsi;

    // object that computes objective and constraints
    FG_eval fg_eval(coeffs);

    //
    // NOTE: You don't have to worry about these options
    //
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
            options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
            constraints_upperbound, fg_eval, solution);

    // Check some of the solution values
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    // TODO: Return the first actuator values. The variables can be accessed with
    // `solution.x[i]`.
    //
    // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
    // creates a 2 element double vector.

    for (int i = 0; i < N - 1; i++) {
        this -> x_predicted[i] = solution.x[x_start + 1 + i];
        this -> y_predicted[i] = solution.x[y_start + 1 + i];
    }

    return {solution.x[delta_start],  solution.x[a_start]};
}
