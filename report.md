## Model Predictive Control Project

This project implements Model Predictive Control (MPC) to drive the car around the track. Unlike in the previous [PID Control project](https://github.com/anguillanneuf/CarND-PID-Control-Project), no cross track error is given this time. Additionally, there's a 100 millisecond latency between actuations commands on top of the connection latency. A Model Predictive Controller is able to drive the car much more smoothly than a PID controller, as can be seen in my project videos for [Smooth Driving at 50 mph](https://youtu.be/6iqupZT016w) and [Smooth Driving at max. 100+ mph](https://youtu.be/9hv8Y8NaNYo). 

[//]: # (Image References)
[image1]: ./img/Screen%20Shot%202017-06-18%20at%203.03.11%20PM.png
[image2]: ./img/Screen%20Shot%202017-06-18%20at%203.03.45%20PM.png
[image3]: ./img/Screen%20Shot%202017-06-18%20at%203.03.52%20PM.png

### The Model

> Student describes their model in detail. This includes the state, actuators and update equations.

The model used for the project is a kinematic model, which are simplifications of dynamic models because it ignores tire forces, gravity, and mass. It involves the followign parameters as the state input [x, y, ψ, ν] and actuator input [δ, a]:  

x is the car's x position in the global coordinate system  
y is the car's y position in the global coordinate system  
ψ (psi) is the car's orientation in radians from the x-direction  
ν is the car's velocity  

δ (delta) is the steering input in radians  
a is the car's accelerator input  

Using these inputs, the MPC solves for the best set of actuator inputs (steering, accelerator) that allow the car to follow a planned trajectory with minimum cost. The cost takes into consideration desired position, heading, stopping, acceleration, steering, change in acceleration and change in steering. Once the car moves, the MPC updates the state at regular intervals and uses it once again to solve future actuator inputs for future trajectories.  

The update equations for state and actuator are in `mpc.cpp` (line 132-168).  

In order to make the car reach very high speeds (80+ mph) on the track, a mechanism must be invoked to help it slow down around the sharp corners. This is done by calculating the [radius of curvature](http://www.intmath.com/applications-differentiation/8-radius-curvature.php), as suggested by Vivek Yadav in the Slack #p-mpc channel.  On top of that, steering angles also inform the MPC what the ideal speed shall be for the next set of actuations. The following chunk of code is inserted into `FG_eval()` to update the ideal driving speed under different conditions. Basically, if the radius of curvature is smaller than 50.0, suggesting a sharp bend, or the steering angles are more than 4 degrees in magnitude, set the ideal speed to half of its original speed.  

```
AD<double> dfdx = coeffs[1] + 2 * coeffs[2] * vars[x_start] + 3 * coeffs[3] * pow(vars[x_start], 2);
AD<double> dfdx2 = 2 * coeffs[2] + 2 * 3 * coeffs[3] * vars[x_start];
AD<double> R_curve = pow(1 + dfdx * dfdx, 1.5) / abs(dfdx2);
double v; 

if (vars[delta_start] > 4.0 * M_PI /180.0 || vars[delta_start] < -4.0 * M_PI /180.0 || R_curve < 50.0){
    v = fmin(55.0, ref_v/2.0);
} else {
    v = fmin(105.0, ref_v);
}
```
MPC applies the brake to slow down the car, turning -7.79˚ travelling at 86.33 mph, around this corner.  

![alt text][image1]

MPC continues to accelerate the car, already travelling at 99+ mph, because both its steering angle and the road curvature are within the safe bounds.  

![alt text][image3]

### Timestep Length and Elapsed Duration (N & dt)

> Student discusses the reasoning behind the chosen N (timestep length) and dt (elapsed duration between timesteps) values. Additionally the student details the previous values tried.  

`N` is the number of timesteps into the future. `dt` is the interval between each activation. Their product gives the predicted trajectory in green. `N` is proportional to the length of this trajectory while `dt` is inversely proportional to it.  

Here are some values that I have tried for `N` (third last) and `dt` (second last) that work for given speeds.  

```
 * "Usage ./mpc cte_cost epsi_cost speed_cost throttle_cost steer_cost dt_throttle_cost dt_steer_cost N dt v"
 * ./mpc 1 1 1 1 500 1 1 10 0.05 50
 * ./mpc 1 1 1 1 500 1 1 20 0.05 55
 * ./mpc 1 1 1 1 100 1 3000 15 0.05 60
 * ./mpc 1 1 1 1 100 1 3000 15 0.05 65
 * ./mpc 1 1 1 0.1 500 5 3000 20 0.05 75
 * ./mpc 1 1e-2 1e-3.5 1e-6 300 1e-2 7000 10 0.05 95
 * ./mpc 1 1 1 1e-6 300 1 7000 12 0.05 100
```
Generally speaking, 0.05 - 0.1 for `dt` and 10 - 20 for `N` are ideal at speeds varying from 50 mph to 70 mph. But `N` > 10 works poorly at higher speeds. The resulted trajectories align with the reference trajectories better farther ahead than immediately ahead, causing the car to adjust itself too aggressively and drive off the track.  

### Polynomial Fitting and MPC Preprocessing

> A polynomial is fitted to waypoints.

> If the student preprocesses waypoints, the vehicle state, and/or actuators prior to the MPC procedure it is described.  

Before fitting the polynomial, the waypoints provided in global coordinate system are first converted to points in the car coordinate system. 

```
// map global x and y positions of waypoints to car coordinate system.
Eigen::VectorXd x(ptsx.size());
Eigen::VectorXd y(ptsy.size());

for (int i = 0; i < ptsx.size(); i++){
    x[i] = (ptsx[i] - px) * cos(psi) + (ptsy[i] - py) * sin(psi);
    y[i] = (ptsy[i] - py) * cos(psi) - (ptsx[i] - px) * sin(psi);
}

// fit a 3rd degree polynomial to the above x and y waypoints
auto coeffs = polyfit(x, y, 3);
```

### Model Predictive Control with Latency

> The student implements Model Predictive Control that handles a 100 millisecond latency. Student provides details on how they deal with latency.

The latency term is applied to each of the state inputs before `cte` and `epsi` are calculated. 

```
// apply latency to state vector
double latency = 0.1;

Eigen::VectorXd state(6);

// initialize to 0.0 because car is at (0,0) in its own coordinate system.
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
```
