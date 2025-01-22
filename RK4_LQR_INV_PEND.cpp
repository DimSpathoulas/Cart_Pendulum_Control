#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// System Parameters
#define M 1.0    // cart mass (kg)
#define m 0.1    // pendulum mass (kg)
#define l 0.5    // pendulum length (m)
#define g 9.81   // gravity (m/s^2)

// Simulation Parameters
#define h 0.01    // time step (s)
#define T 10.0   // total time (s)
#define N 1000   // number of steps

// Control constraints
#define U_MAX 20.0 
#define x_max 8.0
#define x_dot_max 4.0
#define theta_dot_max M_PI/2.0

// State structure
typedef struct {
    double x;      // cart position
    double theta;  // pendulum angle
    double x_dot;  // cart velocity
    double theta_dot;  // pendulum angular velocity
} State;


void compute_accelerations(const State *x_current, double u, double *x_ddot, double *theta_ddot) {
    double sin_theta = sin(x_current->theta);
    double cos_theta = cos(x_current->theta);
    double den = M + m - m * cos_theta * cos_theta;

    // Cart acceleration
    *x_ddot = (u + m * l * x_current->theta_dot * x_current->theta_dot * sin_theta 
               - m * g * sin_theta * cos_theta) / den;

    // pendulum acceleration
    *theta_ddot = (-u * cos_theta - m * l * x_current->theta_dot * x_current->theta_dot * sin_theta * cos_theta 
                   + (M + m) * g * sin_theta) / (l * den);
}


void dynamics_rk4(State *x_next, const State *x_current, double u) {
    State k1, k2, k3, k4;
    double x_ddot, theta_ddot;

    // Enforce control limits on u
    u = fmax(-U_MAX, fmin(U_MAX, u));

    // Compute k1
    compute_accelerations(x_current, u, &x_ddot, &theta_ddot);
    k1.x = x_current->x_dot;
    k1.theta = x_current->theta_dot;
    k1.x_dot = x_ddot;
    k1.theta_dot = theta_ddot;

    // Compute k2
    State temp = *x_current;
    temp.x += h * 0.5 * k1.x;
    temp.theta += h * 0.5 * k1.theta;
    temp.x_dot += h * 0.5 * k1.x_dot;
    temp.theta_dot += h * 0.5 * k1.theta_dot;
    compute_accelerations(&temp, u, &x_ddot, &theta_ddot);
    k2.x = temp.x_dot;
    k2.theta = temp.theta_dot;
    k2.x_dot = x_ddot;
    k2.theta_dot = theta_ddot;

    // Compute k3
    temp = *x_current;
    temp.x += h * 0.5 * k2.x;
    temp.theta += h * 0.5 * k2.theta;
    temp.x_dot += h * 0.5 * k2.x_dot;
    temp.theta_dot += h * 0.5 * k2.theta_dot;
    compute_accelerations(&temp, u, &x_ddot, &theta_ddot);
    k3.x = temp.x_dot;
    k3.theta = temp.theta_dot;
    k3.x_dot = x_ddot;
    k3.theta_dot = theta_ddot;

    // Compute k4
    temp = *x_current;
    temp.x += h * k3.x;
    temp.theta += h * k3.theta;
    temp.x_dot += h * k3.x_dot;
    temp.theta_dot += h * k3.theta_dot;
    compute_accelerations(&temp, u, &x_ddot, &theta_ddot);
    k4.x = temp.x_dot;
    k4.theta = temp.theta_dot;
    k4.x_dot = x_ddot;
    k4.theta_dot = theta_ddot;

    // Final state update
    x_next->x = x_current->x + h / 6.0 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
    x_next->theta = x_current->theta + h / 6.0 * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta);
    x_next->x_dot = x_current->x_dot + h / 6.0 * (k1.x_dot + 2 * k2.x_dot + 2 * k3.x_dot + k4.x_dot);
    x_next->theta_dot = x_current->theta_dot + h / 6.0 * (k1.theta_dot + 2 * k2.theta_dot + 2 * k3.theta_dot + k4.theta_dot);
    
    // Apply state constraints
    x_next->x = fmax(-x_max, fmin(x_next->x, x_max));  // Position limit
    x_next->x_dot = fmax(-x_dot_max, fmin(x_next->x_dot, x_dot_max)); // Velocity limit            
    x_next->theta = atan2(sin(x_next->theta), cos(x_next->theta));  // Normalize angle
    x_next->theta_dot = fmax(-theta_dot_max, fmin(x_next->theta_dot, theta_dot_max)); // Angular velocity limit
}

             
const double K[4] = {1., 30.75724192, 1.97285529, 6.54360696};  // Pre-computed LQR gains
// const double K[4] = {-10.0, -62.63178311, -9.9032684, -13.32610631};
double calc_input_u(const State *current, const State *desired) {
    // Compute state errors
    double x_error = current->x - desired->x;

    double theta_error = current->theta - desired->theta;
    // Normalize angle error to [-pi, pi]
    theta_error = atan2(sin(theta_error), cos(theta_error));

    // double theta_error = atan2(sin(current->theta - M_PI), cos(desired->theta - M_PI));
    
    double x_dot_error = current->x_dot - desired->x_dot;
    double theta_dot_error = current->theta_dot - desired->theta_dot;
    
    // Apply LQR control law: u = -K * (x - x_desired)
    double u = -(K[0] * x_error + 
                 K[1] * theta_error + 
                 K[2] * x_dot_error + 
                 K[3] * theta_dot_error);
    
    // Saturate control input
    return fmax(-U_MAX, fmin(U_MAX, u));
    // return 0.0;
}


int main() {
    // Allocate state arrays
    State *x = (State*)malloc((N+1) * sizeof(State));
    double *u = (double*)malloc(N * sizeof(double));
    
    // Initialize states
    x[0].x = 0;
    x[0].theta = M_PI - 0.95;  // initial deviation
    x[0].x_dot = 0;
    x[0].theta_dot = 0;
    
    // Desired final state
    State x_desired = {0, M_PI, 0, 0};  
    
    printf("%s\n", "Control of a simple pendulum using RK4 and LQR");

    // Simulation loop
    FILE *fp = fopen("trajectory.csv", "w");
    fprintf(fp, "time,x,theta,x_dot,theta_dot,u\n");
    
    for (int k = 0; k < N; k++) {
        // Current time
        double t = k * h;
        
        // Calculate control input
        u[k] = calc_input_u(&x[k], &x_desired);
        
        // Compute next state
        dynamics_rk4(&x[k+1], &x[k], u[k]);
        
        // Save data
        fprintf(fp, "%f,%f,%f,%f,%f,%f\n",
               t,x[k].x, x[k].theta , x[k].x_dot, x[k].theta_dot, u[k]);
    }
    
    fclose(fp);
    free(x);
    free(u);
    return 0;
}
