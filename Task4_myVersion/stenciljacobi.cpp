#define _USE_MATH_DEFINES
// #include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
// #include <regex>
#include <vector>
#include <cmath>
#include <ctime>
#include <assert.h>
#include <omp.h>

using namespace std;

double up_func(double x, double y, double k){
    return sin(k*x)*sinh(k*y);
}

double f(double x, double y, double k){
    return pow(k, 2)*up_func(x, y, k);
}

double getX(size_t i, size_t N, double h){
    return double(i%(N-2)+1)*h;
}

double getY(size_t i, size_t N, double h){
    return (floor(i/(N-2))+1)*h;
}

int jacobiMethod(vector<double>& xk, const vector<double>& b, const double aii, const double aij, const int N, const double maxIter, const int chunk, double& maxTime_s){
    const size_t vec_size = xk.size(), innerLen = N-2;
    size_t iterationsDone = 0;
    bool stop = false;
    double daii = 1/aii, runtime = 0.0, t0 = omp_get_wtime();
    vector<double> xkp1(vec_size, 0); // Vector initialized with 0
    for (size_t iteration = 0; iteration < maxIter; iteration++)
    {
        #pragma omp parallel for schedule(dynamic, chunk) shared(b,xk,xkp1) firstprivate(innerLen,aij,daii,vec_size)
        for (size_t i = 0; i < vec_size; i++) // i is index of vector
        {
            double temp = 0.0;
            if (i>innerLen-1)
            {
                temp = temp - xk[i-innerLen];
            }
            if (i>0 && i%innerLen!=0)
            {
                temp = temp - xk[i-1]; //left of diag
            }
            if (i<(vec_size-1) && (i+1)%innerLen!=0)
            {
                temp = temp - xk[i+1]; //right of diag
            }
            if (i<vec_size-(innerLen))
            {
                temp = temp - xk[i+innerLen];
            }
            xkp1[i] = (temp*aij + b[i])*daii;
        }
        xk.swap(xkp1);
        runtime = omp_get_wtime() - t0;
        stop = runtime > maxTime_s;
        if (stop)
        {
            iterationsDone = iteration+1;
            break;
        }
    }
    if (iterationsDone == 0){iterationsDone = maxIter;}
    maxTime_s = runtime;
    return iterationsDone;
}   
int checkMaxThreads(int& threads){
    if (threads > omp_get_max_threads()){
        threads = omp_get_max_threads();
        return 1;
    }
    return 0;
}
int inputIUE(int argc, char* argv[],size_t &N,int &threads){
    // ./stenciljacobi [resolution] [threads]
    assert(argc==3);
    cout << "Command line arguments recognized. Parsing in IUE mode..." << endl;
    string str = argv[1];
    N = stoi(str);
    str = argv[2];
    threads = stoi(str);
    return 0;
}
int inputPC(int argc, char* argv[], size_t& N, size_t& maxIterations, int& threads)
{   // ./stenciljacobi [resolution] [maxIterations] [threads]
    assert(argc==4);
    cout << "Command line arguments recognized. Parsing in PC mode..." << endl;
    string str = argv[1];
    N = stoi(str);
    str = argv[2];
    maxIterations = stoi(str);
    str = argv[3];
    threads = stoi(str);
    return 1;
}
/*----------------- MAIN -------------------*/
int main(int argc, char *argv[]){
    size_t N, innerLen, vec_size, maxIterations = 200000;
    int threads=0, mode, chunk; 
    const double k= 2 * M_PI, k2 = pow(k,2);
    double h, h2, dh2, aii, aij, maxTime = 10.0;
    //----------------Input----------------//
    // mode = inputPC(argc,argv, N, maxIterations,threads);
    mode = inputIUE(argc,argv,N,threads);
    assert(mode==0 or mode==1);
    cout << "Resolution of <"<< N << "> selected!" << endl;
    cout << "Solving with Jacobi method using <" << maxIterations << "> iterations." << endl; 
    if(checkMaxThreads(threads)){
        cout << "Too many threads selected! Using maximum number instead!" << endl;
    }
    else
    {
        cout << omp_get_max_threads() << " would be available!" << endl;
    }
    cout << "Using <" << threads << "> threads." << endl;
    cout << "Policy: " << "dynamic" << endl;
    // Starting time measurement here
    cout << "-------------------------------------------------------------------" << endl;
    // Let's compute some values we will need quite often
    innerLen = N - 2;
    h = 1.0/double(N-1);
    h2 = pow(h,2);
    dh2 = 1.0/h2;
    vec_size = innerLen*innerLen;
    omp_set_num_threads(threads);
    chunk = int(double(vec_size)/threads);
    //----------------Create FD stencil----------------//
    //Reserve memory
    vector<double> u(vec_size, 0.0), up(vec_size), b(vec_size);
    //Create --> put into function
    #pragma omp parallel for shared(up,b) firstprivate(N,h,k,dh2) schedule(dynamic,chunk)
    for (size_t i = 0; i < vec_size; i++)
    {
        double x = getX(i, N, h);
        double y = getY(i, N, h);
        // Create vectors
        up[i] = up_func(x, y, k);
        if (vec_size-i<=(innerLen)) { // top most row has boundary condition
            b[i] = f(x,y,k) + dh2*up_func(x,1,k); 
        }
        else {
            b[i] = f(x,y,k);
        }
    }

    // Use Jacobi Method to solve for u
    aii = 4.0/h2 + k2;
    aij = -1.0/h2;
    size_t iterationsDone = jacobiMethod(u, b, aii, aij, N, maxIterations, chunk, maxTime);

    cout << "Total # of terations:          " << iterationsDone << endl;
    cout << "Total runtime:                 " << maxTime << "s" << endl;
    cout << "Average runtime per iteration: " << maxTime/iterationsDone << "s/it" << endl;
    return 0;
}
