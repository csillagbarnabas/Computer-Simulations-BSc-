#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <initializer_list>
#include <cmath>
#include <ostream>
#include <sstream>
#include <fstream>

const double pi = 4 * atan(1.0);

template<typename A, typename B, typename C>
std::vector<A> f(std::vector<A>& y,std::vector<A>& ym, const B c,const C cv){
    int n = static_cast<int>(y.size());
    std::vector<A> f(n);
    f[0] = 0;
    f[n-1] = 0;
    for(int i = 1; i < n-1; i++){
        f[i] = 2*y[i] - ym[i] + c*(y[i+1] + y[i-1] - 2*y[i]) / cv;
    }
    return f;
} //egy lépés

template<typename A, typename B, typename C, typename D>
void stepper(std::ofstream& ofile, std::vector<A>& y, B c, C dx, D dt, int m){
    auto cv = dx / dt;
    int n = static_cast<int>(y.size());
    std::vector<A> ym = y; //kezdetben álló húr
    for(int i = 0; i < n; i++){
        ofile << ym[i] << '\t';
    }
    ofile << '\n';
    for(int j = 1; j < m; j++){
        for(int i = 0; i < n; i++){
            ofile << y[i] << '\t';
        }
        ofile << '\n';
        std::vector<A> yp = f(y,ym,c,cv);
        ym = y;
        y = yp;
    }
} //léptető függvény

template<typename T, typename D, typename E>
std::vector<T> kihuzas(const int n, const T A,const D dx,const E o){
    std::vector<T> kihuzas(n+1);
    for(int i = 0; i < int(n/o); i++){
        kihuzas[i] = i*dx*(A/(int(n/o)*dx));
    }
    for(int i = 0; i < static_cast<int>(kihuzas.size()) - int(n/o); i++){
        kihuzas[i + int(n/o)] = A - i*dx*(A/((static_cast<int>(kihuzas.size())-1 - int(n/o))*dx));
    }
    kihuzas[n] = 0;
    return kihuzas;
}//kezdőfeltétel: kihúzott húr

template<typename T, typename D>
std::vector<T> normal(const int n, const T A,const D dx, const int mod){
    std::vector<T> normal(n);
    for(int i = 0; i < n; i++){
        normal[i] = A*sin(i*mod*pi/(n-1));
    }
    return normal;
}

int main(int, char**) {
    double dx = 0.02; //térbeli felosztás
    double dt = 0.02; //időbeli felosztás
    double c = 1; //a hullámok sebessége
    int oa = 3; // a húr hosszát ennyivel osztjuk, és az osztott hossz lesz a fal és a maximum amplitudó távolsága
    int ob = 2;
    double A = 0.13; //maximum amplitudó
    int n = 2501;//beosztások száma a húron +/- 1 KF-től függően
    int mod1 = 5;
    int mod2 = 6;
    double a = 5;
    double b = 8;
    int m = 500;//ennyi időlépést tesz a program
//    std::vector<double> y = kihuzas(n,A,dx,o);
    std::vector<double> ya = normal(n,A,dx,mod1);
    std::vector<double> yb = normal(n,A,dx,mod2);
    std::vector<double> y(static_cast<int>(ya.size()));
    for(int i = 0; i < static_cast<int>(ya.size()); i++){
        y[i] = a*ya[i] + b*yb[i];
    }
//    std::ofstream ofile("kftest3.txt");
//    if(ofile.is_open()){
//        for(int i = 0; i < static_cast<int>(y.size()); i++){
//            ofile << y[i] << '\t';
//        }
//    ofile.close();
//    }
    std::ofstream ofile2("f5b3.dat");
    if(ofile2.is_open()){
        stepper(ofile2,y,c,dx,dt,m);
    ofile2.close();
    }
}