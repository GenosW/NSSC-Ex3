#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <cstring>
#include <cmath>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/src/SparseExtra/MarketIO.h>
#include <omp.h>

using namespace std;


void displayVector(vector<double> x)
{
	for(size_t i = 0; i<x.size(); i++)
	{
		cout << x[i] << endl;
	}
	cout << endl;
}

void displayVectorBig(vector<double> x, string name)
{
	int i = 0;
	while(i < x.size())
	{
		cout << name << "[" << i << "] = " << x[i] << endl;
		i += 100;
	}
	cout << endl;
}


Eigen::SparseMatrix<double> readMatrixEigen(string file)
{
	// this reads in a spd matrix in coordinate form and gives back a full Sparse Matrix

	ifstream in(file);

	string a, b, c, line, entry;
	typedef Eigen::Triplet<double> T;
	vector<T> tripletList;
	int z = 0;
	bool first = true;

	while(getline(in, line))
	{
		if(line.rfind("%",0) != 0)
		{
			stringstream ss(line);
			ss >> a >> b >> c;
			tripletList.push_back(T(stoi(a)-1,stoi(b)-1,stod(c)));
			if(first)
			{
				z = stoi(a);
				first = false;
			}
		}	 
	}
	in.close();

	tripletList.erase(tripletList.begin());

	typedef Eigen::SparseMatrix<double> M;
	M Matrix(z,z);
	Matrix.setFromTriplets(tripletList.begin(),tripletList.end());
	Matrix = M(Matrix.selfadjointView<Eigen::Lower>());

	return Matrix;
}


vector<vector<double>> readMatrix(string file)
{
	// This reads the matrix and converts it to CCS in one run

	ifstream in(file);

	string a, b, c, line, entry;
	vector<vector<double>> result;
	vector<double> V, IA, JA;

	double v, off = 1.0;
	if (file == "test.mtx")
	{
		off = 0.0;
	}

	while(getline(in, line))
	{
		if(line.rfind("%",0) != 0)
		{
			stringstream ss(line);
			ss >> a >> b >> c;
			IA.push_back(stod(a)-off);
			JA.push_back(stod(b)-off);
			V.push_back(stod(c));
		}	 
	}
	in.close();

	// erase first line --> see mtx file
	V.erase(V.begin());
	IA.erase(IA.begin());
	JA.erase(JA.begin());

	// This is where the conversion to CCS starts

	vector<double> JAnew;
	JAnew.push_back(0);
	int i = 0;
	while (i < JA.size())
	{
		double c = 0;
		int currentValue = JA[i];
			while (i+c < JA.size() && currentValue == JA[i + c])
			{
				c++;
			}
		i += c;
		JAnew.push_back(JAnew[JAnew.size() - 1] + c);
	}

	result.push_back(IA);
	result.push_back(JAnew);
	result.push_back(V);
    cout << endl << "Read in everything successfully for our own implementation" << endl << endl;
	return result;
}


void matrixmulCCS(vector<double>& V, vector<double>& IA, vector<double>& JA, vector<double>& x, vector<double>& y)
{
	size_t vec_size = x.size();
	for(size_t i = 0; i<vec_size; i++) 
	{
		for(size_t j = 0; j < i; j++) 
		{
			for(size_t index = JA[j]; index < JA[j+1]; index++)
			{
				if(i == IA[index])
				{
					y[i] += V[index]*x[j];
					break;
				}
			}
		}
		for(size_t index = JA[i]; index < JA[i+1]; index++)
		{
			y[i] += V[index]*x[IA[index]];
		}
	}
}


double scalarProduct(vector<double> a, vector<double> b)
{
	double sum = 0.0;
	for(int i=0; i<a.size(); i++)
	{
		sum += a[i]*b[i];
	}
	return sum;
}
vector<double> scalarVector(double s, vector<double>& v)
{
	vector<double> result;
	for(int i = 0; i<v.size(); i++)
	{
		result.push_back(v[i]*s);
	}

	return result;
}

vector<double> vectorSum(vector<double> a, vector<double> b)
{
	vector<double> result(a.size(),0.0);
	for(int i = 0; i<a.size(); i++)
	{
		result[i] = a[i]+b[i];
	}

	return result;
}

vector<double> vectorDiff(vector<double> a, vector<double> b)
{
	vector<double> result(a.size(),0.0);
	for(int i = 0; i<a.size(); i++)
	{
		result[i] = a[i]-b[i];
	}

	return result;
}

double euclidicNorm(vector<double> a)
{
	return sqrt(scalarProduct(a,a));
}


vector<vector<double>> CG(vector<vector<double>> A, vector<double> x0, vector<double> r0, int iteration, vector<double> xorig)
{
	ofstream outfile;
	outfile.open("data.txt");
	outfile << "||rk||2 / ||r0||2" << endl;
	int vec_size = A[1].size() - 1;
	vector<vector<double>> result;
	vector<double> IA, JA, V;
	IA = A[0];
	JA = A[1];
	V = A[2];
	vector<double> r = r0;
	vector<double> p = r0;
	vector<double> x = x0;
	double beta;
	double rkrk;
	double alpha;

	for(int i = 0; i<iteration; i++)
	{
		vector<double> Ap(JA.size(), 0.0);
		matrixmulCCS(V,IA,JA,p,Ap);
		rkrk = scalarProduct(r,r);
		alpha = rkrk/(scalarProduct(p,Ap));
		x = vectorSum(x,scalarVector(alpha,p));
		r = vectorDiff(r,scalarVector(alpha,Ap));
		beta = scalarProduct(r,r)/rkrk;
		p = vectorSum(r,scalarVector(beta,p));
		if(i%20 == 0)
		{
			cout << "Own implementation iteration " << i << " done" << endl << endl ;
			outfile << euclidicNorm(r)/euclidicNorm(r0) << endl;
		}
		
	}

	result.push_back(x);
	result.push_back(r);
	outfile.close();
	return result;
}

int main(int argc, char *argv[])
{
	string filename = argv[1]; 
	vector<vector<double>> A = readMatrix(filename);
	int vec_size = A[1].size() - 1;
    
	vector<double> xstar(vec_size, 1.0);
	vector<double> b(vec_size, 0.0);		// right hand side for x = [1 1 1 ... 1]

	matrixmulCCS(A[2],A[0],A[1],xstar,b);

	vector<double> x0(vec_size, 0.0);
	vector<vector<double>> res = CG(A,x0,b,stoi(argv[2]),xstar);
	vector<double> xnew = res[0];

	cout << "After "<< argv[2]  << " iterations, solution vector x equals: " << endl << endl;
	displayVectorBig(xnew, "x");
	cout << "||rk||2 / ||r0|||2 = " << euclidicNorm(res[1])/euclidicNorm(b) << endl << endl;
    


    // Eigen starts here

    Eigen::setNbThreads(2);
    Eigen::SparseMatrix<double> SparseMat = readMatrixEigen(filename);
	Eigen::VectorXd xStarE, xE(vec_size);
	xStarE.setOnes(vec_size);
	Eigen::VectorXd bE = SparseMat * xStarE;
	

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper > cg(SparseMat);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double> > cgCholesky(SparseMat);

    ofstream outfileEigen;
    ofstream outfileEigenCholesky;
	outfileEigen.open("dataEigen.txt");
    outfileEigen << "Error for diagonally preconditioned Eigen method" << endl;
    outfileEigenCholesky.open("dataEigenCholesky.txt");
	outfileEigenCholesky << "Error for Cholesky preconditioned Eigen method" << endl;

    
    for(int i = 0; i<stoi(argv[2]); i+=20)
    {
        cg.setMaxIterations(i);
		cg.setTolerance(1e-15);
        xE = cg.solve(bE);
        if(i%100 == 0)
        {
            cout << "Diagonally preconditioned Eigen iteration " << i << endl;
        }
        outfileEigen << cg.error() << endl;
    }
	outfileEigen.close();
	
    
    for(int i = 0; i<stoi(argv[2]); i+=20)
    {
        cgCholesky.setMaxIterations(i);
		cgCholesky.setTolerance(1e-15);
        xE = cgCholesky.solve(bE);
        if(i%100 == 0)
        {
            cout << "Incomplete Cholesky preconditioned Eigen iteration " << i << endl;
        }
        outfileEigenCholesky << cgCholesky.error() << endl;
    }
    outfileEigenCholesky.close();

	cout << endl << endl << endl;
	cout << "------------ Summary after " << argv[2] << " iterations ------------" << endl << endl;
	cout << "||rk||2 / ||r0|||2 for our own implementation = " << euclidicNorm(res[1])/euclidicNorm(b) << endl << endl;
	cout << "Error of diagonally preconditioned Eigen method =" << cg.error() << endl << endl;
	cout << "Error of Incomplete Cholesky preconditioned Eigen method =" << cgCholesky.error() << endl << endl;;
    

    return 0;
}