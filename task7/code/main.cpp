#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

template<typename A>
std::vector<A> stepper(std::vector<A> const & v, int N, int n){//élő határ
    std::vector<A> stepper = v;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            int m = 0;
            if (i+1<N && v[i+1+j*N] == 1){
                m += 1;
            }
            if (j+1 < N && v[i+(j+1)*N] == 1){
                m += 1;
            }
            if (i+1<N && j+1 < N  && v[i+1+(j+1)*N] == 1){
                m += 1;
            }
            if (i-1 > -1  && v[i-1+j*N] == 1){
                m += 1;
            }
            if (j-1 > -1  &&v[i+(j-1)*N] == 1){
                m += 1;
            }
            if (j-1 > -1  && i+1<N && v[i+1+(j-1)*N] == 1){
                m += 1;
            }
            if (i-1 > -1  && j+1<N &&v[i-1+(j+1)*N] == 1){
                m += 1;
            }
            if (i-1 > -1  && j-1>-1 &&v[i-1+(j-1)*N] == 1){
                m += 1;
            }
            if (m==n+1){
                stepper[i+N*j] = 1;
            }
            if (m<n){
                stepper[i+N*j] = 0;
            }
            if (m>n+1){
                stepper[i+N*j] = 0;
            }
        }
    }
    return stepper;
}

template<typename A>
std::vector<A> stepper2(std::vector<A> const & v, int N, int n){//konstans véletlen határ
    std::vector<A> stepper2 = v;
    for(int i = 1; i < N-1; i++){
        for (int j = 1; j < N-1; j++){
            int m = 0;
            if (i+1<N-1 && v[i+1+j*N] == 1){
                m += 1;
            }
            if (j+1 < N-1 && v[i+(j+1)*N] == 1){
                m += 1;
            }
            if (i+1<N-1 && j+1 < N-1  && v[i+1+(j+1)*N] == 1){
                m += 1;
            }
            if (i-1 > 0  && v[i-1+j*N] == 1){
                m += 1;
            }
            if (j-1 > 0  &&v[i+(j-1)*N] == 1){
                m += 1;
            }
            if (j-1 > 0  && i+1<N-1 && v[i+1+(j-1)*N] == 1){
                m += 1;
            }
            if (i-1 > 0  && j+1<N-1 &&v[i-1+(j+1)*N] == 1){
                m += 1;
            }
            if (i-1 > 0  && j-1>0 &&v[i-1+(j-1)*N] == 1){
                m += 1;
            }
            if (m==n+1){
                stepper2[i+N*j] = 1;
            }
            if (m<n){
                stepper2[i+N*j] = 0;
            }
            if (m>n+1){
                stepper2[i+N*j] = 0;
            }
        }
    }
    return stepper2;
}

template<typename A>
std::vector<A> stepper3(std::vector<A> const & v, int N, int n){//periodikus határfeltétel
    std::vector<A> stepper3 = v;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            int m = 0;

            if (i+1<N && v[i+1+j*N] == 1){
                m += 1;
            }
            if (i+1==N && v[j*N] == 1){
                m += 1;
            }

            if (j+1 < N && v[i+(j+1)*N] == 1){
                m += 1;
            }
            if (j+1 == N && v[i] == 1){
                m += 1;
            }

            if (i+1<N && j+1 < N  && v[i+1+(j+1)*N] == 1){
                m += 1;
            }
            if (i+1 == N && j+1 == N  && v[0] == 1){
                m += 1;
            }

            if (i-1 > -1  && v[i-1+j*N] == 1){
                m += 1;
            }
            if (i-1 == -1  && v[N-1+j*N] == 1){
                m += 1;
            }

            if (j-1 > -1  &&v[i+(j-1)*N] == 1){
                m += 1;
            }
            if (j-1 == -1  &&v[i+(N-1)*N] == 1){
                m += 1;
            }

            if (j-1 > -1  && i+1<N && v[i+1+(j-1)*N] == 1){
                m += 1;
            }
            if (j-1 == -1  && i+1==N && v[(N-1)*N] == 1){
                m += 1;
            }

            if (i-1 > -1  && j+1<N &&v[i-1+(j+1)*N] == 1){
                m += 1;
            }
            if (i-1 == -1  && j+1==N &&v[N-1] == 1){
                m += 1;
            }

            if (i-1 > -1  && j-1>-1 &&v[i-1+(j-1)*N] == 1){
                m += 1;
            }
            if (i-1 == -1  && j-1==-1 &&v[N-1+(N-1)*N] == 1){
                m += 1;
            }
            if (m==n+1){
                stepper3[i+N*j] = 1;
            }
            if (m<n){
                stepper3[i+N*j] = 0;
            }
            if (m>n+1){
                stepper3[i+N*j] = 0;
            }
        }
    }
    return stepper3;
}

template<typename A>
void initializer(std::vector<A> & v){
    for(int i = 0; i < v.size(); i++){
        v[i] = rand() % 2;
//        std::cout << v[i] << '\t';
    }
}

int main(int, char**) {
    int N = 200;
    std::vector<int> T(N*N);
    int M = 400;
    int n = 2;
    initializer(T);
    std::string out = "conwaypern2h.dat";
    std::ofstream dataFile(out);
    for(int k = 0; k<M;k++){
        for(int i = 0; i < T.size(); i++){
            dataFile << T[i]<< '\t';
        }
        dataFile << '\n';
        T = stepper3(T,N,n);
    }
}