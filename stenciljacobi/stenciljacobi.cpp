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

void displayVector(const vector<double>& v, size_t N){
    int limit = 16;
    if (limit > N){limit = N;}
    for (size_t i = 0; i < limit; i++)
    {
        cout << "\n" << setw(4) << v[i];
    }
    cout << endl;
}

void displayMatrix(vector<double>& A, size_t N){
    int limit = 10;
    if (limit > N){limit = N;}
    cout << "mat:" << setw(4) << "\n";
    for (size_t i = 0; i < limit; i++)
    {
        for (size_t j = 0; j < limit; j++)
        {
            cout << A[i*N + j] << "\t";
        }
        cout << endl;
    }
}

void displayGrid(vector<double>& v, size_t N){
    int limit = 10;
    if (limit > N){limit = N;}
    cout << setw(8) << setprecision(4) << "\n" ;
    for (size_t i = 0; i < limit; i++)
    {
        for (size_t j = 0; j < limit; j++)
        {
            cout  << v[i*N + j] << "\t\t";
        }
        cout << "\n";
    }
}

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

int vecVecSub(vector<double>& c,vector<double>& a, vector<double>& b, size_t vec_size){
    assert(vec_size == c.size() && vec_size == a.size() && vec_size == b.size());
    for (size_t i = 0; i < vec_size; i++)
    {
        c[i]=(a[i] - b[i]);
    }
    return 0;
}

double normVecMax(const vector<double>& v, size_t vec_size, size_t& index){
    assert(vec_size == v.size());
    double temp = 0;
    index = -1;
    for (size_t i = 0; i < vec_size; i++)
    {
        if (abs(v[i])>temp)
        {
            temp = abs(v[i]);
            index = i;
        }
    }
    return temp;
}

double normVecEuc(const vector<double>& v, size_t vec_size){
    assert(vec_size == v.size());
    double temp = 0;
    for (size_t i = 0; i < vec_size; i++)
    {
        temp += pow(v[i],2);
    }
    // return sqrt(temp/v_size); //"normalized variant"
    return sqrt(temp);
}
int initVec(vector<double>& u, size_t vec_size){
    assert(vec_size == u.size());
    for (size_t i = 0; i < vec_size; i++)
    {
        u[i] = 0;
    }
    return 0;
}

int checkEq(vector<double>& result, vector<double>& u, vector<double>& b, double aii, double aij, int N){
    const size_t vec_size=u.size(), innerLen = N-2;
    if (vec_size!=pow(innerLen,2)){
        cerr << "vec_size != pow(N,2)" << endl;
        return -1;
    }
    for (size_t i = 0; i < vec_size; i++)
    {
        double tempDiff = aii*u[i] - b[i];
        if (i>(innerLen-1)){
            tempDiff += aij*u[i-innerLen];
        }
        if ((i>0) && (i%innerLen!=0)){
            tempDiff += aij*u[i-1];
        }
        if (i<(vec_size-1) && (i+1)%innerLen!=0){
            tempDiff += aij*u[i+1];
        }
        if (i<vec_size-innerLen){
            tempDiff += aij*u[i+innerLen];
        }
        result[i]= tempDiff; 
    }
    return 0;
}

int jacobiMethod(vector<double>& xk, const vector<double>& b, const double aii, const double aij, const int N, const double maxIter){
    const size_t vec_size = xk.size(), innerLen = N-2;
    double temp=0.0;
    int chunk = innerLen;
    vector<double> xkp1(vec_size, 0); // Vector initialized with 0
    //# pragma omp for schedule(dynamic, chunk) nowait shared(aij,aii,N,h,k,dh2,up,b)
    for (size_t iteration = 0; iteration < maxIter; iteration++)
    {
        // xkp1 = (b[i] - aij*xk[j])/aii
        // # pragma omp sections nowait
        // # pragma omp section
        for (size_t i = 0; i < vec_size; i++) // i is index of vector
        {
            temp = b[i];
            if (i>innerLen-1)
            {
                temp = temp - aij*xk[i-innerLen];
            }
            if (i>0 && i%innerLen!=0)
            {
                temp = temp - aij*xk[i-1]; //left of diag
            }
            if (i<(vec_size-1) && (i+1)%innerLen!=0)
            {
                temp = temp - aij*xk[i+1]; //right of diag
            }
            if (i<vec_size-(innerLen))
            {
                temp = temp - aij*xk[i+innerLen];
            }
            xkp1[i] = temp/aii;
        }
        // # pragma omp section
        xk = xkp1;
        // for (size_t i = 0; i < vec_size; i++)
        // {
        //     xk[i] = xkp1[i];
        // }
    }
    return 0;
}     
/*----------------- MAIN -------------------*/
int main(int argc, char *argv[]){
    size_t N, innerLen, vec_size, maxIterations;
    const double k= 2 * M_PI, k2 = pow(k,2);
    double h, h2, dh2, eucNorm, maxNorm, aii, aij, itRuntime = 0.0;
	struct timespec tItEnd, tItStart;
    //----------------Input----------------//
    assert(argc==3);
    cout << "Command line arguments recognized. Parsing..." << endl;
    string str = argv[1];
    N = stoi(str);
    str = argv[2];
    maxIterations = stoi(str);
    cout << "Resolution of <"<< N << "> selected!" << endl;
    cout << "Solving with Jacobi method using <" << maxIterations << "> iterations." << endl; 
    
    // Starting time measurement here
    cout << "-------------------------------------------------------------------" << endl;
    // Let's compute some values we will need quite often
    innerLen = N - 2;
    h = 1.0/double(N-1);
    h2 = pow(h,2);
    dh2 = 1.0/h2;
    //----------------Create FD stencil----------------//
    //Reserve memory
    vec_size = pow(innerLen,2);
    //
    vector<double> u(vec_size, 0.0), up(vec_size, 0.0), b(vec_size, 0.0);
    omp_set_num_threads(2);
    int chunk = N-2;
    //Create --> put into function
    #pragma omp for schedule(dynamic, chunk) nowait shared(N,h,k,dh2,up,b)
    for (size_t i = 0; i < vec_size; i++)
    {
        // Compute x and y depending on index i
        // x = double(i%(innerLen)+1)*h;
        // y = (floor(i/(innerLen))+1)*h;
        double x = getX(i, N, h);
        double y = getY(i, N, h);
        // Create vectors
        up[i] = up_func(x, y, k);
        if (vec_size-i<=(innerLen)) { // top most row has boundary condition
            b[i] = f(x, y, k) + dh2*up_func(x,1,k); 
        }
        else {
            b[i] = f(x, y, k);
        }
    }
    
    // Use Jacobi Method to solve for u
    aii = 4.0/h2 + k2;
    aij = -1.0/h2;
    clock_gettime(CLOCK_REALTIME, &tItStart);
    jacobiMethod(u, b, aii, aij, N, maxIterations);
    clock_gettime(CLOCK_REALTIME, &tItEnd);

    vector<double> result(vec_size, 0.0);
    // Compute residual = A*u - b after Jacobi Method
    checkEq(result, u, b, aii, aij, N);
    size_t index_of_max=0;
    eucNorm = normVecEuc(result, vec_size);
    maxNorm = normVecMax(result, vec_size, index_of_max);
    cout << "eucNorm(A*u-b): " << eucNorm << endl;
    cout << "maxNorm(A*u-b): " << maxNorm << endl;

    // Compute  total error up - u 
    vecVecSub(result, up, u, vec_size); // result = up - u
    eucNorm = normVecEuc(result, vec_size);
    maxNorm = normVecMax(result, vec_size, index_of_max);
    cout << "eucNorm(up-u): " << eucNorm << endl;
    cout << "maxNorm(up-u): " << maxNorm << endl;

    // Get runtime for this run and write into file
    itRuntime = double(tItEnd.tv_sec - tItStart.tv_sec) + double(tItEnd.tv_nsec - tItStart.tv_nsec)/1e9;
    cout << "runtime: " << itRuntime << "s" << endl;
    cout << "Iterations: " << maxIterations << endl;
    return 0;
}

/*
resolution={16, 32, 64, 128, 256}
*/