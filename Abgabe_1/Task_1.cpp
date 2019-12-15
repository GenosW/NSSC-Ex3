#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <cstring>
#include <cmath>

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
    cout << endl << "Read in everything successfully" << endl << endl;
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
	outfile << "||rk||2 / ||r0|||2;||ek||A" << endl;
	int vec_size = A[1].size() - 1;
	vector<vector<double>> result;
	vector<double> IA, JA, V;
	IA = A[0];
	JA = A[1];
	V = A[2];
	vector<double> r = r0;		// starting paramter
	vector<double> p = r0;		// starting paramter
	vector<double> x = x0;		// starting paramter
	double beta;
	double rkrk;
	double alpha;

	for(int i = 0; i<iteration; i++)
	{
		vector<double> Ap(JA.size(), 0.0);
		matrixmulCCS(V,IA,JA,p,Ap);
		rkrk = scalarProduct(r,r);
		alpha = rkrk/(scalarProduct(p,Ap));			// step length
		x = vectorSum(x,scalarVector(alpha,p));		// new approximate solution
		r = vectorDiff(r,scalarVector(alpha,Ap));	// residual
		beta = scalarProduct(r,r)/rkrk;				// improvement this step
		p = vectorSum(r,scalarVector(beta,p));		// search direction for next iteration
		if(i%10 == 0)
		{
			cout << "Iteration " << i << " done" << endl << endl ;
		}
		// This is for outputting the two norms. if only end result is needed, one can comment the next four lines out
		vector<double> ek = vectorDiff(xorig, x);
		vector<double> Aek(vec_size, 0.0);
		matrixmulCCS(A[2],A[0], A[1], ek, Aek);
		outfile << euclidicNorm(r)/euclidicNorm(r0) << ";" << sqrt(scalarProduct(ek,Aek)) << endl;
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
	vector<double> x(vec_size, 1.0);
	vector<double> b(vec_size, 0.0);		// right hand side for x = [1 1 1 ... 1]

	matrixmulCCS(A[2],A[0],A[1],x,b);

	vector<double> x0(vec_size, 0.0);
	vector<vector<double>> res = CG(A,x0,b,stoi(argv[2]),x);
	vector<double> xnew = res[0];

	cout << "After "<< argv[2]  << " iterations, solution vector x equals: " << endl << endl;
	displayVectorBig(xnew, "x");
	cout << "||rk||2 / ||r0|||2 = " << euclidicNorm(res[1])/euclidicNorm(b) << endl << endl;
	vector<double> ek = vectorDiff(x, xnew);
	vector<double> Aek(vec_size, 0.0);
	matrixmulCCS(A[2],A[0], A[1], ek, Aek);
	cout << "||ek||A = " << sqrt(scalarProduct(ek,Aek)) << endl;

    return 0;
}
