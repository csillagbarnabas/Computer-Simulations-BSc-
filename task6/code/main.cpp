#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip> 

#include "RK.h"        // ODE integration routines, Runge-Kutta ...

const double pi = 4 * atan(1.0);
const double r1 = 1.1;
const double r2 = 1.1;

const double alpha = 2.2;
const double beta = 2.1;
const double k1 = 0.9;
const double k2 = 1.5;

const double a = 0.8;
const double b = 0.002;
const double c = 0.001;
const double d = 0.8;

std::vector<double> logi(const std::vector<double>& y){
    double t = y[0];
    double x = y[1];
    std::vector<double> logi={1,1.1*x*(1-x)};
    return logi;
}

std::vector<double> csatlog1(const std::vector<double>& y1,const std::vector<double>& y2){
    double t = y1[0];
    double n1 = y1[1];
    double n2 = y2[1];
    std::vector<double> csatlog1 = {1,r1*n1*(1-(n1+alpha*n2)/k1)};
    return csatlog1;
}

std::vector<double> csatlog2(const std::vector<double>& y1,const std::vector<double>& y2){
    double t = y1[0];
    double n1 = y1[1];
    double n2 = y2[1];
    std::vector<double> csatlog2 = {1,r2*n2*(1-(n2+beta*n1)/k2)};
    return csatlog2;
}

std::vector<double> LV1(const std::vector<double>& y1,const std::vector<double>& y2){
    double t = y1[0];
    double n1 = y1[1];
    double n2 = y2[1];
    std::vector<double> LV1 = {1,a*n1-b*n1*n2};
    return LV1;
}

std::vector<double> LV2(const std::vector<double>& y1,const std::vector<double>& y2){
    double t = y1[0];
    double n1 = y1[1];
    double n2 = y2[1];
    std::vector<double> LV2 = {1,c*n1*n2-d*n2};
    return LV2;
}

int main() {
    double t = 0;
//    double x = 1.0;
    double n1 = 500;
    double n2 = 500;
    std::vector<double> y1 = {t,n1};
    std::vector<double> y2 = {t,n2};
    double tMax = 100;
    double dt = 0.0005;
    double accuracy = 1e-6;
    std::string out = "lotka8.dat";
    std::ofstream dataFile(out);
    while (t < tMax) {

//        RK4Step(x, dt, f5);
//        t = x[0], theta = x[1], omega = x[2];
//        Euler_fo(y,dt, logi);
//        adaptiveRK4Step(y, dt,accuracy, logi);
        Euler_fo_2csat(y1,y2,dt, LV1,LV2);
//        RK4Step_2(y1,y2,dt, LV1, LV2);
        t = y1[0], n1 = y1[1];
        n2 = y2[1];
        dataFile << std::setprecision(16) << t << '\t' << n1 << '\t' << n2 << '\n';
    }
    std::cout << n1 << '\t' << n2<< '\n';
//    std::cout << "[ " << k1 <<"," << '\t' << k2 <<"," << '\t' <<alpha <<"," << '\t' << beta <<"]," << '\n';
    std::cout << " Output data to file " + out << std::endl;
    dataFile.close();
}