#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
}


int main(int argc, const char *argv[]) {
    uWS::Hub h;

    // speed up testing by entering weights for each term in the cost function
    /*
     * ./mpc 1 1 1 1 500 1 1 10 0.05 50
     * ./mpc 1 1 1 1 500 1 1 10 0.05 55
     * ./mpc 1 1 1 1 100 1 3000 15 0.05 65
     * ./mpc 1 1 1 1 100 1 3000 15 0.05 75
     * ./mpc 1 1 1 0.1 500 5 3000 20 0.05 75
     * ./mpc 1 1e-2 1e-3.5 1e-6 300 1e-2 7000 10 0.05 100
     * ./mpc 1 1 1 1e-6 300 1 7000 12 0.05 85
     * */

    double w0, w1, w2, w3, w4, w5, w6, dt_input, v_input;
    int N_input;

    if (argc == 11 ) {
        w0 = strtod( argv[1], NULL ) ;
        w1 = strtod( argv[2], NULL ) ;
        w2 = strtod( argv[3], NULL ) ;
        w3 = strtod( argv[4], NULL ) ;
        w4 = strtod( argv[5], NULL ) ;
        w5 = strtod( argv[6], NULL ) ;
        w6 = strtod( argv[7], NULL ) ;
        N_input = strtod( argv[8], NULL ) ;
        dt_input = strtod( argv[9], NULL ) ;
        v_input = strtod( argv[10], NULL ) ;
    } else {
        cout << "Usage ./mpc w0 w1 w2 w3 w4 w5 w6 N dt ref_v" << endl;
        return -1 ;
    }

    // MPC is initialized here!
    MPC mpc;
    mpc.Init(w0, w1, w2, w3, w4, w5, w6, N_input, dt_input, v_input);


    h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                       uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        string sdata = string(data).substr(0, length);
        cout << sdata << endl;
        if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
            string s = hasData(sdata);
            if (s != "") {
                auto j = json::parse(s);
                string event = j[0].get<string>();
                if (event == "telemetry") {
                    // j[1] is the data JSON object
                    vector<double> ptsx = j[1]["ptsx"];
                    vector<double> ptsy = j[1]["ptsy"];
                    double px = j[1]["x"];
                    double py = j[1]["y"];
                    double psi = j[1]["psi"];
                    double v = j[1]["speed"];
                    double steer = j[1]["steering_angle"]; // in radians
                    double throttle = j[1]["throttle"];

                    /*
                    * TODO: Calculate steering angle and throttle using MPC.
                    *
                    * Both are in between [-1, 1].
                    *
                    */

                    // Map global x and y positions of waypoints to car coordinate system.
                    Eigen::VectorXd x(ptsx.size());
                    Eigen::VectorXd y(ptsy.size());

                    for (int i = 0; i < ptsx.size(); i++){
                        x[i] = (ptsx[i] - px) * cos(psi) + (ptsy[i] - py) * sin(psi);
                        y[i] = (ptsy[i] - py) * cos(psi) - (ptsx[i] - px) * sin(psi);
                    }

                    // fit a 3rd degree polynomial to the above x and y waypoints
                    auto coeffs = polyfit(x, y, 3);

                    // apply latency to state vector
                    double latency = 0.1;

                    Eigen::VectorXd state(6);

                    // initialize to 0.0 because we are in car coordinate system.
                    px = py = psi = 0.0;

                    // update px, py, psi, v according to latency
                    px += v * cos(psi) * latency;
                    py += v * sin(psi) * latency;
                    psi -= v * steer / mpc.getLf() * latency;
                    v += throttle * latency;

                    // calculate cross track error;
                    double cte = polyeval(coeffs, px) - py;

                    // calculate orientation error;
                    double epsi = psi - atan(coeffs[1] + 2 * coeffs[2] * px + 3 * coeffs[3] * pow(px, 2));

                    state << px, py, psi, v, cte, epsi;

                    // solve for steer and throttle given state and coeffs
                    double steer_value;
                    double throttle_value;

                    auto vars = mpc.Solve(state, coeffs);

                    // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
                    // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
                    steer_value = -vars[0]/deg2rad(25);
                    throttle_value = vars[1];

                    steer_value = fmin(fmax(-1.0, steer_value), 1.0);
                    throttle_value = fmin(fmax(-1.0, throttle_value), 1.0);

                    json msgJson;

                    msgJson["steering_angle"] = steer_value;
                    msgJson["throttle"] = throttle_value;

                    //Display the MPC predicted trajectory
                    vector<double> mpc_x_vals;
                    vector<double> mpc_y_vals;

                    //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                    // the points in the simulator are connected by a Green line
                    mpc_x_vals = mpc.x_predicted;
                    mpc_y_vals = mpc.y_predicted;

                    msgJson["mpc_x"] = mpc_x_vals;
                    msgJson["mpc_y"] = mpc_y_vals;

                    //Display the waypoints/reference line
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;

                    //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                    // the points in the simulator are connected by a Yellow line

                    /*for (int i = 1; i < x.size(); i++)
                    {
                        next_x_vals.push_back(x[i]);
                        next_y_vals.push_back(y[i]);
                    };*/

                    for(int i = 0; i < 100; i += 5) {
                        next_x_vals.push_back(i * 1.0);
                        next_y_vals.push_back(polyeval(coeffs, i));
                    }

                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;


                    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
                    std::cout << msg << std::endl;
                    // Latency
                    // The purpose is to mimic real driving conditions where
                    // the car does actuate the commands instantly.
                    //
                    // Feel free to play around with this value but should be to drive
                    // around the track with 100ms latency.
                    //
                    // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
                    // SUBMITTING.
                    this_thread::sleep_for(chrono::milliseconds(100));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
