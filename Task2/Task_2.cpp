#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <cstring>
#include <assert.h>
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
}
vector<vector<double>> readMatrix(string file)
{
	ifstream in(file);

	string a, b, c, line, entry;
	vector<vector<double>> result;
	vector<double> V, IA, JA;
	double off = 1.0;
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
	// erase first line --> see mtx file
	V.erase(V.begin());
	IA.erase(IA.begin());
	JA.erase(JA.begin());

	vector<double> JAnew;
	JAnew.push_back(0);
	uint i = 0;
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
	cout << "size(IA): " << result[0].size() << endl;
	//displayVector(IA);
	result.push_back(JAnew);
	cout << "size(JA): " << result[1].size()  << endl;
	//displayVector(JAnew);
	result.push_back(V);
	cout << "size(V): " << result[2].size()  << endl;
	//displayVector(V);
    cout << "Read in everything successfully" << endl;
	return result;
}
void v_matrixmulCCS(vector<double>& V, vector<double>& IA, vector<double>& JA, vector<double>& x, vector<double>& y)
{
	size_t vec_size = x.size();
	//double sum;
	for(size_t i = 0; i<vec_size; i++) // iterate over all rows
	{
		for(size_t j = 0; j < i; j++) // iterate over all columns left of diag
		{
			for(size_t index = JA[j]; index < JA[j+1]; index++)
			{
				if(i == IA[index])
				{
					cout << "y[i] :" << y[i] << endl;
					y[i] += V[index]*x[j];
					cout << "i: " << i << " index: " << index << " V[index]: " << V[index] << " y[i]: " << y[i] << endl;
					break;
				}
			}
		}
		for(size_t index = JA[i]; index < JA[i+1]; index++)
		{
			cout << "y[i] :" << y[i] << endl;
			y[i] += V[index]*x[IA[index]];
			// y.at(i) += V.at(index)*x.at(IA.at(index));
			cout << "i: " << i << " index: " << index << " V[index]: " << V[index] << " x[index]: " << x[IA[index]] << " y[i]: " << y[i] << endl;
		}
		cout << "y[" << i << "]= " << y[i] << endl;
	}
}
void matrixmulCCS(vector<double>& V, vector<double>& IA, vector<double>& JA, vector<double>& x, vector<double>& y)
{
	size_t vec_size = x.size();
	//double sum;
	for(size_t i = 0; i<vec_size; i++) // iterate over all rows
	{
		for(size_t j = 0; j < i; j++) // iterate over all columns left of diag
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
		//if(i%1000 == 0) cout << "y[" << i << "]= " << y[i] << endl;
	}
}
int input(int argc, char* argv[], string& mode,string& filename,int& threads){
    assert(argc==4);
	// ./cgEigen [resolution] [filename] [threads]
    cout << "Command line arguments recognized. Parsing..." << endl;
    mode = argv[1];
	cout << "Mode: <" << mode << "> selected." << endl;
	if (argc > 1){
		filename = argv[2];
		cout << "Loading file <" << filename << ">..." << endl;
		if (argc > 2){
			threads = stoi(argv[3]);
			cout << "Using <" << threads << "> threads." << endl;
		}
	}
	return 0;
}

int main(int argc, char *argv[])
{
	string filename = "test.mtx"; // "s3rmt3m1.mtx"
	string mode = "std";
	int threads = 4;
	if (argc != 0)
	{
		input(argc,argv,mode,filename,threads);
	}

	Eigen::setNbThreads(threads);

	vector<vector<double>> A = readMatrix(filename);
	int vec_size = A[1].size() - 1;
	vector<double> xStar(vec_size, 1.0);
	vector<double> b(vec_size, 0.0);
	if (filename == "test.mtx"){
		v_matrixmulCCS(A[2],A[0],A[1],xStar,b);
	}
	else
	{
		matrixmulCCS(A[2],A[0],A[1],xStar,b);
	}

	Eigen::SparseMatrix<double> SpA;//(vec_size+1,vec_size+1);
	cout << "here" << endl;
	Eigen::loadMarket(SpA,filename);
	cout << "here" << endl;
	Eigen::VectorXd xStarE, bE(vec_size), xE(vec_size), x0;
	xStarE.setOnes(vec_size);
	x0.setZero(vec_size);

	//cout << xStarE << endl;
	// fill A and b
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> cg(SpA);
	bE = SpA * xStarE;
	// cg.setMaxIterations(1);
	//for(int i=0; i < 1000; i++)
	// {
	// 	cg.setMaxIterations(i);
	// 	xE = cg.solve(bE);
	// 	if (i%100 == 0)
	// 	{
	// 		cout << "iterations this loop iteration: " << cg.iterations() << endl;
	// 		cout << cg.error() << endl;;
	// 		cout << (xE - xStarE).norm() << endl;;
	// 		cout << "-------" << endl;
	// 	}
	// }
	cg.setMaxIterations(100000);
	xE = cg.solveWithGuess(bE,x0);
	cout << "-------" << endl;
	std::cout << "#iterations:			" << cg.iterations() << std::endl;
	std::cout << "estimated error:		" << cg.error() << std::endl;
	std::cout << "Norm(x - x*):			" << (xE - xStarE).norm() << std::endl;
	std::cout << "NormA(x-x*):			" << (xE-xStarE).transpose()*(SpA*(xE-xStarE)) << std::endl;
	std::cout << "maxCoeffs(x - x*):    " << (xE - xStarE).maxCoeff() << std::endl;
	// update b, and solve again
	//x = cg.solve(b);
// DiagonalPreconditioner
// IncompleteCholesky
    return 0;
}