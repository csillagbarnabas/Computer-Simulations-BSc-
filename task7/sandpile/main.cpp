#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

template<typename A>
std::vector<A> stepper_per(std::vector<A> const & v, int N){//periodikus határfeltétel
    std::vector<A> steppe_perr = v;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (v[i + j*N]>3){
                if(i < N-1 && j < N-1 && i >0 && j >0){
                    stepper_per[i+1+j*N] += 1;
                    stepper_per[i-1+j*N] += 1;
                    stepper_per[i+(j+1)*N] += 1;
                    stepper_per[i+(j-1)*N] += 1;
                    stepper_per[i+j*N] += -4;
                }
                if(i == N-1 && j < N-1 && i >0 && j >0){
                    stepper_per[j*N] += 1;
                    stepper_per[i-1+j*N] += 1;
                    stepper_per[i+(j+1)*N] += 1;
                    stepper_per[i+(j-1)*N] += 1;
                    stepper_per[i+j*N] += -4;
                }
                if(i < N-1 && j == N-1 && i >0 && j >0){
                    stepper_per[j*N] += 1;
                    stepper_per[i-1+j*N] += 1;
                    stepper_per[i] += 1;
                    stepper_per[i+(j-1)*N] += 1;
                    stepper_per[i+j*N] += -4;
                }
                if(i < N-1 && j < N-1 && i == 0 && j >0){
                    stepper_per[i+1+j*N] += 1;
                    stepper_per[N-1+j*N] += 1;
                    stepper_per[i+(j+1)*N] += 1;
                    stepper_per[i+(j-1)*N] += 1;
                    stepper_per[i+j*N] += -4;
                }
                if(i < N-1 && j < N-1 && i >0 && j == 0){
                    stepper_per[i+1+j*N] += 1;
                    stepper_per[i-1+j*N] += 1;
                    stepper_per[i+(j+1)*N] += 1;
                    stepper_per[i+(N-1)*N] += 1;
                    stepper_per[i+j*N] += -4;
                }
            }
        }
    }
    return stepper_per;
}

template<typename A>
std::vector<A> stepper_open(std::vector<A> const & v, int N){//nyílt határfeltétel
    std::vector<A> stepper_open = v;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (v[i + j*N]>3){
                if(i < N-1 && j < N-1 && i >0 && j >0){
                    stepper_open[i+1+j*N] += 1;
                    stepper_open[i-1+j*N] += 1;
                    stepper_open[i+(j+1)*N] += 1;
                    stepper_open[i+(j-1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i == N-1 && j < N-1 && i >0 && j >0){
                    stepper_open[i-1+j*N] += 1;
                    stepper_open[i+(j+1)*N] += 1;
                    stepper_open[i+(j-1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i < N-1 && j == N-1 && i >0 && j >0){
                    stepper_open[j*N] += 1;
                    stepper_open[i-1+j*N] += 1;
                    stepper_open[i+(j-1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i < N-1 && j < N-1 && i == 0 && j >0){
                    stepper_open[i+1+j*N] += 1;
                    stepper_open[i+(j+1)*N] += 1;
                    stepper_open[i+(j-1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i < N-1 && j < N-1 && i >0 && j == 0){
                    stepper_open[i+1+j*N] += 1;
                    stepper_open[i-1+j*N] += 1;
                    stepper_open[i+(j+1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i < N-1 && j < N-1 && i == 0 && j == 0){
                    stepper_open[i+1+j*N] += 1;
                    stepper_open[i+(j+1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i == N-1 && j < N-1 && i > 0 && j == 0){
                    stepper_open[i-1+j*N] += 1;
                    stepper_open[i+(j+1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i < N-1 && j == N-1 && i == 0 && j >0){
                    stepper_open[i+1+j*N] += 1;
                    stepper_open[i+(j-1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
                if(i == N-1 && j == N-1 && i >0 && j >0){
                    stepper_open[i-1+j*N] += 1;
                    stepper_open[i+(j-1)*N] += 1;
                    stepper_open[i+j*N] += -4;
                }
            }
        }
    }
    return stepper_open;
}

template<typename A>
void initializer(std::vector<A> & v){
    for(int i = 0; i < v.size(); i++){
        v[i] = 7;
    }
}

int main(int, char**) {
    int N = 100;
    std::vector<int> T(N*N);
    int M = 10000;
    int i1 = 50;
    int j1 = 50;
    initializer(T);
    std::string out = "sandpileflow.dat";
    std::ofstream dataFile(out);
    for(int k = 0; k<M;k++){
        if(k % 50 == 0){
        for(int i = 0; i < T.size(); i++){
            dataFile << T[i]<< '\t';
        }
        dataFile << '\n';
        }
        T = stepper_open(T,N);
        T[i1 + j1*N] += 3; //plusz ráöntés egy pontban
    }
}
