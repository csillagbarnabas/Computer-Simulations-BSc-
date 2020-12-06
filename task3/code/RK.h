#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

template<typename T, typename F>
std::vector<T> multiply_vs(std::vector<T> & x, const F y){
    std::vector<T> res(x.size());
    for(int i = 0; i<x.size();i++){
        res[i] = x[i]*y;
    }
    return res;
}
template<typename T, typename F>
std::vector<T> add(std::vector<T> & x,std::vector<F> & y){
    std::vector<T> res(x.size());
    for(int i = 0; i<x.size();i++){
        res[i] = x[i]+y[i];
    }
    return res;
}

template<typename T, typename F,typename RHS>
void EulerCromer (std::vector<T>& x,const F dt, RHS f) {
    x[2] = x[2] + f(x)[2]*dt;
    x[1] = x[1] + x[2]*dt;
    x[0] = x[0] + dt;
}

template<typename T, typename F,typename RHS>
void Euler(std::vector<T>& x,const F dt, RHS f) {
    std::vector<T> a = f(x);
//    x[2] = x[2] + f(x)[2]*dt;
//    x[1] = x[1] + f(x)[1]*dt; //sajnos csak így működik
//    x[0] = x[0] + dt;
     x = {x[0] + dt,x[1] + dt*a[1],x[2]+dt*a[2]};
}

//  Fourth order Runge-Kutta
template<typename T, typename F, typename RHS>
void RK4Step(std::vector<T>& x, F tau, RHS derivs)
{
//    std::vector<T> k1 = tau * derivs(x);
    std::vector<T> perm = derivs(x);
    std::vector<T> k1 = multiply_vs(derivs(x),tau);
    std::vector<T> k2 = multiply_vs(derivs(add(x,multiply_vs(k1,0.5))),tau);
    std::vector<T> k3 = multiply_vs(derivs(add(x,multiply_vs(k2,0.5))),tau);
    std::vector<T> k4 = multiply_vs(derivs(add(x,k3)),tau);
    x = add(x,multiply_vs(add(add(k1, multiply_vs(k2,2)), add(multiply_vs(k3,2), k4)),1.0/6.0));
//    x = {x[0] + (k1[0]+k2[0]*2+k3[0]*2+k4[0])/6,x[1] + (k1[1]+k2[1]*2+k3[1]*2+k4[1])/6,x[2] + (k1[2]+k2[2]*2+k3[2]*2+k4[2])/6};
}
//  Adaptive step size control using Runge-Kutta and step doubling
template<typename T, typename F, typename L, typename RHS>
void adaptiveRK4Step(std::vector<T>& x, F tau, L accuracy, RHS derivs)
{
    const double SAFETY = 0.9, PGROW = -0.2, PSHRINK = -0.25,
                 ERRCON = 1.89E-4, TINY = 1.0e-30;
    int n = int(x.size());
    std::vector<T> x_half(n), x_full(n), Delta(n);
    std::vector<T> scale = derivs(x);
    for (int i = 0; i < n; i++)
        scale[i] = abs(x[i]) + abs(scale[i] * tau) + TINY;
    double err_max;
    while (true) {
        // take two half steps
        double tau_half = tau / 2;
        x_half = x;
        RK4Step(x_half, tau_half, derivs);
        RK4Step(x_half, tau_half, derivs);
        // take full step
        x_full = x;
        RK4Step(x_full, tau, derivs);
        // estimate error
        std::vector<double> Delta(n);
        for(int i = 0; i < n; i++){
            Delta[i] = x_half[i] - x_full[i];
        }
        err_max = 0;
        for (int i = 0; i < n; i++)
            err_max = std::max(err_max, std::abs(Delta[i]) / scale[i]);
        err_max /= accuracy;
        if (err_max <= 1.0)
            break;
        double tau_temp = SAFETY * tau * pow(err_max, PSHRINK);
        if (tau >= 0.0)
            tau = std::max(tau_temp, 0.1 * tau);
        else
            tau = std::min(tau_temp, 0.1 * tau);
        if (std::abs(tau) == 0.0) {
            std::cerr << "adaptiveRK4Step: step size underflow\naborting ..."
                 << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    tau *= (err_max > ERRCON ? SAFETY * pow(err_max, PGROW) : 5.0);
    for(int i=0; i < n; i++){
        x[i] = x_half[i] + Delta[i] / 15.0;
    }
}