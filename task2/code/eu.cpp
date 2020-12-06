#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "vector.hpp"        
#include "odeint.hpp"        
#include "odeint.cpp"
#include "vector.cpp"

const double pi = 4 * atan(1.0);

const double g = 9.8;        // acceleration of gravity

const double L = 12;              // length of pendulum
const double q = 1;              // damping coefficient
const double Omega_D = 1;    // frequency of driving force
const double F_D = 3;            // amplitude of driving force
const double m = 5;             //single pendulum weight

const double L1 = 1; //double pendulum lengths
const double L2 = 1;
const double m1 = 1;
const double m2 = 1;

double theta1 = 0.5236;
double theta2 = theta1 /2;
double omega1 = 0.3;
double omega2 = 0.6;

const std::string single_or_double = "single";

const std::string model = "physical"; 
            /*physical model desc:
                simple - math pendulum
                damped - damped pend
                damped-driven - damped-driven pend
                physical - physical pend
             */

const std::string intg_method = "Euler";
                        /* integration method:
                            Euler
                            Euler-Cromer
                            RK - simple Runge-Kutta 
                            RKA - Runge-Kutta /w adaptive step size
                         */


//Hyperparamters of integration

const double tMax = 20;     
double dt = 0.05;
double accuracy = 1e-6;


double t = 0;
double omega = 4;
double theta = 2;




cpl::Vector f(const cpl::Vector& x) // extended derivative vector
{  
    t = x[0];
    theta = x[1];
    omega = x[2];
    cpl::Vector f(3);             // Vector with 3 components
    f[0] = 1;
    f[1] = omega;
    if (model == "physical")
        f[2] = - (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t);
    if (model == "damped-driven")
        f[2] = - (g/L) * theta - q * omega + F_D * sin(Omega_D * t);
    if (model == "damped")
        f[2] = -(g/L) * theta - q *omega;
    if (model == "simple")
        f[2] = -(g/L) * theta;
    return f;
}

double kin_energy()
{
    return 0.5 * L * L * omega * omega;  
}

double pot_energy()
{
    return (L-L * cos(theta)) * g;
}

void RK45A(cpl::Vector& x)
{
    std::ofstream dataFile("RK45A/" + model + ".dat");

    while (t < tMax)
    {
        
        adaptiveRK4Step(x, dt, accuracy, f);
        t = x[0], theta = x[1], omega = x[2];
        if (model=="physical")
        {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }
        dataFile << t << '\t' << theta << '\t' << omega << '\t' << kin_energy() <<'\t' << pot_energy()<< '\n';
    }
    dataFile.close();
}

void RK4(cpl::Vector& x)
{
    std::ofstream dataFile("RK4/" + model + ".dat");

    while (t < tMax)
    {
        RK4Step(x, dt, f);
        t = x[0], theta = x[1], omega = x[2];
        if (model=="physical")
        {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }
        dataFile << t << '\t' << theta << '\t' << omega << '\t' << kin_energy() <<'\t' << pot_energy()<< '\n';
    }
    dataFile.close();
}

void Euler(cpl::Vector& x)
{
    std::ofstream dataFile("Euler/" + model + ".dat");
    
    while(t < tMax)
    {
        x += dt * f(x);
        t = x[0], theta = x[1], omega = x[2];

        if (model=="physical")
        {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }
        dataFile << t << '\t' << theta << '\t' << omega << '\t' << kin_energy() <<'\t' << pot_energy()<< '\n';
    }
    dataFile.close();
}

void Euler_Cromer(cpl::Vector& x)
{
    std::ofstream dataFile("Euler-Cromer/" + model + ".dat");
    cpl::Vector f(3);
    while(t < tMax)
    {
        
        x[0] += dt;
        t = x[0];
        if (model == "physical")
            x[2] += dt*(- (g/L) * sin(theta) - q * omega + F_D * sin(Omega_D * t));
        if (model == "damped-driven")
            x[2] += dt*(- (g/L) * theta - q * omega + F_D * sin(Omega_D * t));
        if (model == "damped")
            x[2] += dt*(-(g/L) * theta - q *omega);
        if (model == "simple")
            x[2] += dt*(-(g/L) * theta);

        omega = x[2];

        x[1] += omega * dt;
        theta = x[1];
        

        if (model=="physical")
        {
            while (theta >= pi) theta -= 2 * pi;
            while (theta < -pi) theta += 2 * pi;
        }
        dataFile << t << '\t' << theta << '\t' << omega << '\t' << kin_energy() <<'\t' << pot_energy()<< '\n';
    }
    dataFile.close();
}

cpl::Vector derivs_2_pend(const cpl::Vector& y)
{
    cpl::Vector f(5);

    f[0] = y[1];
    double ad = y[2] - y[0]; //difference of angles
    double denominator1 = (m1 + m2) * L1 - m2 * L1 * cos(ad) * cos(ad);
    f[1] = (m2 * L1 * y[1] * y[1] * sin(ad) * cos (ad) +
            m2 * g * sin(y[2]) * cos(ad) +
            m2 * L2 * y[3] * y[3] *sin(ad) -
            (m1 + m2) * g * sin(y[0]) ) / denominator1;
    
    f[2] = y[4];
    
    double denominator2 = (L2/L1) * denominator1;
    f[3] = (-m2 * L2 * y[3] *y[3] * sin(ad) * cos(ad) + 
            (m1 + m2) * g * sin(y[0]) * cos(ad) -
            (m1 + m2) * L1 * y[1] * y[1] * sin(ad) -
            (m1 + m2) * g * sin(y[2])) / denominator2;

    f[4] = 1;
    return f;
}

void int_dp (cpl::Vector& y)
{
    std::ofstream dataFile("double_pendulum/dp.dat");

    while(y[4] < tMax)  
    {
        adaptiveRK4Step(y, dt, accuracy,derivs_2_pend);
        dataFile << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] 
                 << "\t" << y[4] << "\n";
    }
    dataFile.close();
}

int main() 
{
    if (single_or_double == "single")
    {
        cpl::Vector x(3);
        x[0] = t;
        x[1] = theta;
        x[2] = omega;

        if (intg_method == "RK45A"){RK45A(x); return 0;}

        if (intg_method == "RK4"){RK4(x); return 0;}

        if (intg_method == "Euler_Cromer"){Euler_Cromer(x); return 0;}

        if (intg_method == "Euler"){Euler(x); return 0;}
    }

    if (single_or_double == "double")
    {
        cpl::Vector y(5);
        y[0] = t;
        y[1] = theta1;
        y[2] = omega1;
        y[3] = theta2;
        y[4] = omega2;

        int_dp(y);
        return 0;
    }
    

    std::cout << " Done." << std::endl;
    
    return 0;    
}