#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>

#include "RK.h"

const double pi = 4 * atan(1.0);
const double GmPlusM = 4 * pi * pi;
const double alpha = 1.1e-11;

bool switch_t_with_y = false;    //  to interpolate to y = 0

//  Derivative vector for Newton's law of gravitation
std::vector<double> f(const std::vector<double>& x) {
    double t = x[0], r_x = x[1], r_y = x[2], v_x = x[3], v_y = x[4];
    double rSquared = r_x*r_x + r_y*r_y;
    double rCubed = rSquared * std::sqrt(rSquared);
    std::vector<double> f(5);
    f[0] = 1;
    f[1] = v_x;
    f[2] = v_y;
    f[3] = - GmPlusM * r_x / rCubed;
    f[4] = - GmPlusM * r_y / rCubed;
    if (switch_t_with_y) {
        //  use y as independent variable
        for (int i = 0; i < 5; i++)
            f[i] /= v_y;
    }
    return f;
}

std::vector<double> g(const std::vector<double>& x) {
    double t = x[0], r_x = x[1], r_y = x[2], v_x = x[3], v_y = x[4];
    double rSquared = r_x*r_x + r_y*r_y;
    double rCubed = rSquared * std::sqrt(rSquared);
    std::vector<double> g(5);
    g[0] = 1;
    g[1] = v_x;
    g[2] = v_y;
    g[3] = - GmPlusM * r_x * (1 + alpha / rCubed) / rCubed;
    g[4] = - GmPlusM * r_y * (1 + alpha / rCubed) / rCubed;
    if (switch_t_with_y) {
        //  use y as independent variable
        for (int i = 0; i < 5; i++)
            g[i] /= v_y;
    }
    return g;
}

//  Change independent variable from t to y and step back to y = 0
void interpolate_crossing(std::vector<double> x, int& crossing) {
    ++crossing;
    switch_t_with_y = true;
    RK4Step(x, -x[2], f);
    std::cout << " crossing " << crossing << "\t t = " << x[0]
         << "\t x = " << x[1] << std::endl;
    switch_t_with_y = false;
}

int main() {
    std::cout << " Kepler orbit comparing fixed and adaptive Runge-Kutta\n" << std::endl;
    double r_ap = 0.46669835; //aphéliumból indítom
    double eccentricity = 0.2056; //excentritása a Merkúr pályájának
    double a = r_ap / (1 + eccentricity);
    double T = pow(a, 1.5);
    double v0 = sqrt(GmPlusM * (2 / r_ap - 1 / a));
    double periods = 10000.0;
    double dt = 0.001;
    double accuracy = 1e-6;
    std::vector<double> x0(5);
    x0[0] = 0;  x0[1] = r_ap;  x0[2] = 0;  x0[3] = 0;  x0[4] = v0;
    std::ofstream dataFile("f3fixed.data");
    std::vector<double> x = x0;
    int steps = 0, crossing = 0;
    std::cout << "\n Integrating with fixed step size" << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
  //  do {
  //      for (int i = 0; i < 5; i++)
  //          dataFile << x[i] << '\t';
  //      dataFile << '\n';
  //      double y = x[2];
  //      RK4Step(x, dt, f);
  //      ++steps;
//        if (y * x[2] < 0)
//            interpolate_crossing(x, crossing);
  //  } while (x[0] < periods * T);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "It took " << double((static_cast<std::chrono::duration<double, std::milli>>(t2-t1)).count()) <<" milliseconds"<<std::endl;
    std::cout << " number of fixed size steps = " << steps << std::endl;
    std::cout << " data in file fixed.data" << std::endl;
    dataFile.close();
//    dt = 0.05;
    dataFile.open("f3adaptive2.data");
    x = x0;
    steps = crossing = 0;
    double dt_max = 0, dt_min = 100;
    std::cout << "\n Integrating with adaptive step size" << std::endl;
    t1 = std::chrono::high_resolution_clock::now();
    do {
        for (int i = 0; i < 5; i++)
            dataFile << x[i] << '\t';
        dataFile << '\n';
        double t_save = x[0];
        double y = x[2];
        adaptiveRK4Step(x, dt, accuracy, g);
        double step_size = x[0] - t_save;
        ++steps;
//        if (step_size < dt_min) dt_min = step_size;
//        if (step_size > dt_max) dt_max = step_size;
//        if (y * x[2] < 0)
//            interpolate_crossing(x, crossing);
    } while (x[0] < periods * T);
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "It took " << double((static_cast<std::chrono::duration<double, std::milli>>(t2-t1)).count()) <<" milliseconds"<<std::endl;
    std::cout << " number of adaptive steps = " << steps << std::endl;
    std::cout << " step size: min = " << dt_min << "  max = " << dt_max << std::endl;
    std::cout << " data in file adaptive.data" << std::endl;
    dataFile.close();
}