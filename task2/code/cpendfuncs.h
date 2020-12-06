#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "..\headers\vector3.h"
#include <fstream>
#include <list>
template<typename T>
auto changedp(Vector3<T> const& r){//x,y,z
    auto R = sqrt(sqv(r.x)+sqv(r.y)+sqv(r.z));
    return Vector3<T>{R,acos(r.z/R),atan2(r.y,r.x)};
}
template<typename T>
auto changepd(Vector3<T> const& r){//r,phi,theta
    auto sy = sin(r.y);
    return Vector3<T>{r.x*sy*cos(r.z),r.x*sy*sin(r.z),r.x*cos(r.y)};//x,y,z
}
template<typename T>
auto dcomponentssumphi(Vector3<T> const& r1, Vector3<T> const& r2){
    return r1.x+r1.y+r1.z-r2.x-r2.y-r2.z;
}
template<typename T>
auto dcomponentssumtheta(Vector3<T> const& r1, Vector3<T> const& r2){
    return r1.x+r1.y-r2.x-r2.y;
}
template<typename T, typename F, typename C>
auto cpendcasumphi(std::vector<T> const& q,std::vector<Vector3<F>> const& rq, Vector3<F> const& r, C const& q0){
    T sum = (T)0;
    for(int i = 0; i < rq.size(); ++i){
        sum += dcomponentssumphi(rq[i],r) * q0 * q[i]/(length(rq[i]-r)*sqlength(rq[i]-r));
    }
    return sum;
}
template<typename T, typename F, typename C> //külön tesztre készült fv, ami csak az új részeket tartalmazza
auto cpendcasum(std::vector<T> const& q,std::vector<Vector3<F>> const& rq, Vector3<F> const& r, C const& q0){
    T sum = (T)0;
    for(int i = 0; i < rq.size(); ++i){
        sum += q0 * q[i]/(length(rq[i]-r)*sqlength(rq[i]-r));
    }
    return sum;
}
template<typename T, typename F, typename C>
auto cpendcasumtheta(std::vector<T> const& q,std::vector<Vector3<F>> const& rq, Vector3<F> const& r, C const& q0){
    T sum = (T)0;
    for(int i = 0; i < rq.size(); ++i){
        sum += dcomponentssumtheta(rq[i],r) * q0 * q[i]/(length(rq[i]-r)*sqlength(rq[i]-r));
    }
    return sum;
}
template<typename State, typename C, typename kg, typename M>
auto cpenda( Vector3<State> const& r, Vector3<State> const& v, std::vector<Vector3<double>> const& rq, C const& q0,
std::vector<C> const& q,kg const& m, M const& l){
    Vector3<State> rd{changepd(r)};
    auto sy = sin(r.y); auto cy = cos(r.y); auto cz = cos(r.z); auto sz = sin(r.z);
    return Vector3<State>{(State)0,//9.81 is the gravitational acceleration in Hungary
    sy*cy*sqv(v.z)+(double)9.81*sy/l - cpendcasumphi(q,rq,rd,q0)*(cy*cz+cy*sz-sy)/(m*l),
    -2*cy*sqv(r.z)/sy-cpendcasumtheta(q,rq,rd,q0)*(sz-cz)/(m*l*sy)};
}
template<typename State, typename T, typename C, typename kg, typename M>
void cpendsolve_rk4_txt(Vector3<State> const& r0,Vector3<State> const& v0, std::vector<Vector3<double>> const& rq,
T const& t0,T const& t1, int const& N, C const& q0,std::vector<C> const& q,kg const& m, M const& l,
std::ofstream& ttxt,std::ofstream& xtxt,std::ofstream& ytxt,std::ofstream& ztxt,
std::ofstream& vxtxt,std::ofstream& vytxt,std::ofstream& vztxt){
    T h = (t1-t0)/N;
    std::vector<Vector3<State>> r(N+1);
    r[0]=changedp(r0);
    std::vector<Vector3<State>> v(N+1);
    v[0]=v0;
    std::vector<T> t(N+1);
    t[0] = t0;
    int i = 0;
    while(i < N){
        Vector3<State> kv1 = cpenda(r[i],v[i],rq,q0,q,m,l); //mivel az egyenletekben nincs t, ezért nem is léptetem
        Vector3<State> kv2 = cpenda(r[i] + (h * (T)0.5) * kv1,v[i] + (h * (T)0.5) * kv1,rq,q0,q,m,l);
        Vector3<State> kv3 = cpenda(r[i] + (h * (T)0.5) * kv2,v[i] + (h * (T)0.5) * kv2,rq,q0,q,m,l);
        Vector3<State> kv4 = cpenda(r[i] + h * kv3,v[i] + h * kv3,rq,q0,q,m,l);
        v[i+1] = v[i] + (kv1 + kv4 + (T)2 * (kv2 + kv3)) * (h / (T)6);
        Vector3<State> kr1 = v[i];
        Vector3<State> kr2 = v[i] + (h * (T)0.5) * kr1;
        Vector3<State> kr3 = v[i] + (h * (T)0.5) * kr2;
        Vector3<State> kr4 = v[i] + h * kr3;
        r[i+1] = r[i] + (kr1 + kr4 + (T)2 * (kr2 + kr3)) * (h / (T)6);
        t[i+1] = t[i] + h;
        ++i;
    }
    for(int j = 0; j < N; ++j){
        xtxt << changepd(r[j]).x <<std::endl;
        ytxt << changepd(r[j]).y <<std::endl;
        ztxt << changepd(r[j]).z <<std::endl;
        vxtxt << v[j].x <<std::endl;
        vytxt << v[j].y <<std::endl;
        vztxt << v[j].z <<std::endl;
        ttxt << t[j] <<std::endl;
    }
}
template<typename State, typename T, typename C, typename kg, typename M>
void cpendsolve_rk4(Vector3<State> const& r0,Vector3<State> const& v0, std::vector<Vector3<double>> const& rq,
T const& t0,T const& t1, int const& N, C const& q0,std::vector<C> const& q,kg const& m, M const& l,
std::vector<Vector3<State>> & r, std::vector<T> & t){
    T h = (t1-t0)/N;
    r[0]=changedp(r0);
    std::vector<Vector3<State>> v(N+1);
    v[0]=v0;
    t[0] = t0;
    int i = 0;
    while(i < N){
        Vector3<State> kv1 = cpenda(r[i],v[i],rq,q0,q,m,l); //mivel az egyenletekben nincs t, ezért nem is léptetem
        Vector3<State> kv2 = cpenda(r[i] + (h * (T)0.5) * kv1,v[i] + (h * (T)0.5) * kv1,rq,q0,q,m,l);
        Vector3<State> kv3 = cpenda(r[i] + (h * (T)0.5) * kv2,v[i] + (h * (T)0.5) * kv2,rq,q0,q,m,l);
        Vector3<State> kv4 = cpenda(r[i] + h * kv3,v[i] + h * kv3,rq,q0,q,m,l);
        v[i+1] = v[i] + (kv1 + kv4 + (T)2 * (kv2 + kv3)) * (h / (T)6);
        Vector3<State> kr1 = v[i];
        Vector3<State> kr2 = v[i] + (h * (T)0.5) * kr1;
        Vector3<State> kr3 = v[i] + (h * (T)0.5) * kr2;
        Vector3<State> kr4 = v[i] + h * kr3;
        r[i+1] = r[i] + (kr1 + kr4 + (T)2 * (kr2 + kr3)) * (h / (T)6);
        t[i+1] = t[i] + h;
        ++i;
    }
    for(int j=0;j<N;++j){
		r[j] = changepd(r[j]);
	}
}
template<typename State, typename T, typename C, typename kg, typename M, typename RHS>
void solve_rk4(State const& x0,State const& vx0,T const& t0,T const& t1, int const& N,
C k,kg m,M gamma,std::ofstream& ttxt,std::ofstream& xtxt, RHS f){
    T h = (t1-t0)/N;
    std::vector<State> x(N+1);
    x[0]=x0;
    std::vector<State> v(N+1);
    v[0]=vx0;
    std::vector<T> t(N+1);
    t[0] = t0;
    int i = 0;
    while(i < N){
        State kv1 = f(x[i],v[i]);
        State kv2 = f(x[i] + (h * (T)0.5) * kv1,v[i] + (h * (T)0.5) * kv1);
        State kv3 = f(x[i] + (h * (T)0.5) * kv2,v[i] + (h * (T)0.5) * kv2);
        State kv4 = f(x[i] + h * kv3,v[i] + h * kv3);
        v[i+1] = v[i] + (kv1 + kv4 + (T)2 * (kv2 + kv3)) * (h / (T)6);
        State kr1 = v[i];
        State kr2 = v[i] + (h * (T)0.5) * kr1;
        State kr3 = v[i] + (h * (T)0.5) * kr2;
        State kr4 = v[i] + h * kr3;
        x[i+1] = x[i] + (kr1 + kr4 + (T)2 * (kr2 + kr3)) * (h / (T)6);
        t[i+1] = t[i] + h;
        ++i;
    }
    for(int j = 0; j < N; ++j){
        xtxt << x[j] <<std::endl;
        ttxt << t[j] <<std::endl;
    }
}