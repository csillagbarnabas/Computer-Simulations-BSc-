#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip> 

#include "RK.h"        // ODE integration routines, Runge-Kutta ...

const double pi = 4 * atan(1.0);

const double g = 9.8;        // acceleration of gravity

double L = 9.8;              // length of pendulum
double q = 0.2;              // damping coefficient
double Omega_D = 2.0/3.0;    // frequency of driving force
double F_D = 1.2;            // amplitude of driving force
bool nonlinear;              // linear if false

double l2 = 9.8;

std::vector<double> f1(const std::vector<double>& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    std::vector<double> f1(3);             // Vector with 3 components
    f1[0] = 1;
    f1[1] = omega;
    f1[2] = - (g/L) * theta;
    return f1;
}

std::vector<double> f2(const std::vector<double>& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    std::vector<double> f2(3);             // Vector with 3 components
    f2[0] = 1;
    f2[1] = omega;
    f2[2] = - (g/L) * theta - q * omega;
    return f2;
}

std::vector<double> f3(const std::vector<double>& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    std::vector<double> f3(3);             // Vector with 3 components
    f3[0] = 1;
    f3[1] = omega;
    f3[2] = - (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    return f3;
}

std::vector<double> f4(const std::vector<double>& x) {  // extended derivative vector
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    std::vector<double> f4(3);             // Vector with 3 components
    f4[0] = 1;
    f4[1] = omega;
    f4[2] = - (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    return f4;
}

std::vector<double> f5(const std::vector<double>& x){
    double t = x[0];
    double theta = x[1];
    double omega = x[2];
    double theta2 = x[3];
    double omega2 = x[4];
    std::vector<double> f5(5);
    f5[0] = 1;
    f5[1] = omega;
    f5[2] = (l2*omega2*omega2*sin(theta-theta2) + L*omega*omega*sin(theta-theta2)*cos(theta-theta2)+2*g*sin(theta)-g*sin(theta2)*cos(theta-theta2))/(L*(2-cos(theta-theta2)*cos(theta-theta2)));
    f5[3] = omega2;
    f5[4] = (l2*omega2*omega2*sin(theta-theta2)*cos(theta-theta2)+2*L*omega*omega*sin(theta-theta2)+2*g*sin(theta)*cos(theta-theta2)+2*g*sin(theta2))/(l2*(2-cos(theta-theta2)*cos(theta-theta2)));
    return f5;
}

int main() {
    std::cout << " Nonlinear damped driven pendulum\n";
    std::string response = "n";
    nonlinear = true;
    double theta = 3;
    double omega = 0.0;
    double theta2 = 2.8;
    double omega2 = 0.0;
    double tMax = 50;
    double dt = 0.02;
    double accuracy = 1e-6;
    std::string out = "2pendtesta28.dat";
    std::ofstream dataFile(out);

    double t = 0;
//    std::vector<double> x(3);
    std::vector<double> x(5);
    x[0] = t;
    x[1] = theta;
    x[2] = omega;
    x[3] = theta2;
    x[4] = omega2;
    while (t < tMax) {
//        Euler(x,dt,f4);
//        EulerCromer (x,dt,f4);
//        RK4Step(x, dt, f5);
        adaptiveRK4Step(x, dt,accuracy, f5);
//        t = x[0], theta = x[1], omega = x[2];
        t = x[0], theta = x[1], omega = x[2]; theta2 = x[3], omega2 = x[4];
        while (theta >= pi) theta -= 2 * pi;
        while (theta < -pi) theta += 2 * pi;
        while (theta2 >= pi) theta2 -= 2 * pi;
        while (theta2 < -pi) theta2 += 2 * pi;
        dataFile << std::setprecision(16) << t << '\t' << theta << '\t' << omega << '\t' << theta2 << '\t' << omega2 << '\n';
    }

    std::cout << " Output data to file " + out << std::endl;
    dataFile.close();
}

