#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
using namespace std;


string fileName;       // name of output file

void getInput();                 // for user to input parameters
void EulerCromer(double dt);     // takes an Euler-Cromer step
double energy();                 // computes the energy

void EulerCromer (double dt, double a, double v, double x, double omega) {
    a = - omega * omega * x;
    v = v + a * dt;
    x += v * dt;
}

double energy (double v, double x, double omega) {
    return 0.5 * (v * v + omega * omega * x * x);
}

double x_anal(double t, double x0,double v0, double w){
    return x0*cos(w*t) + v0*sin(w*t)/w;
}

double v_anal(double t, double x0,double v0, double w){
    return -x0*w*sin(w*t) + v0*cos(w*t);
}

void anal(double pi, double stepsPerPeriod, double t, int periods, double dt, double a, double v, double x, double omega){
    ofstream ofile("otthonteszte2_anal.txt");
    if (!ofile) {
        cerr << "Cannot open " << fileName << "\nExiting ...\n";
    }
    double T = 2 * pi / omega;
    dt = T / stepsPerPeriod;
    t = 0;
    const double x0 = 30;
    const double v0 = 40;
    for (int p = 1; p <= periods; p++) {
        for (int s = 0; s < stepsPerPeriod; s++) {
            x = x_anal(t,x0,v0,omega);
            v = v_anal(t,x0,v0,omega);
            ofile << t << '\t' << x << '\t' << v << '\n';
            t += dt;
        }
        cout << "Period = " << p << "\tt = " << t
             << "\tx = " << x << "\tv = " << v;
//             << "\tenergy = " << energy() << endl;
    }
    ofile.close();
}
int main ( ) {
    double omega = 13;          // the natural frequency
    double x =30;
    double v = 40;
    double a;           // position and velocity at time t
    int periods = 10;           // number of periods to integrate
    const double pi = 4 * atan(1.0);


//    int stepsPerPeriod = 100;
//    ofstream file("otthonteszte3.txt");
//    if (!file) {
//        cerr << "Cannot open " << fileName << "\nExiting ...\n";
//        return 1;
//    }
    ofstream ofile("idomeres.txt");
    if (!ofile) {
        cerr << "Cannot open " << fileName << "\nExiting ...\n";
    }
    for (int i = 1; i < 100; i++){
    
        double T = 2 * pi / omega;
        double dt = T / (i*1000000);
        double t = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
//    file << t << '\t' << x << '\t' << v << '\n';
        for (int p = 1; p <= 1; p++) {
            for (int s = 0; s < (i*1000000); s++) {
                t += dt;
                a = - omega * omega * x;
                v += a * dt;
                x += v * dt;
//            file << t << '\t' << x << '\t' << v << '\n';
        }
//        cout << "Period = " << p << "\tt = " << t
//             << "\tx = " << x << "\tv = " << v
//             << "\tenergy = " << energy() << endl;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        ofile << double((static_cast<std::chrono::duration<double, std::milli>>(t2-t1)).count()) << '\n';
    }
    ofile.close();
//    file.close();
    cout << "Done"<< endl;
}