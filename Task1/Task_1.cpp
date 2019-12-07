#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <cstring>

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
	// erase first line --> see mtx file
	V.erase(V.begin());
	IA.erase(IA.begin());
	JA.erase(JA.begin());

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
	cout << "size(IA): " << result[0].size() << endl;
	displayVector(IA);
	result.push_back(JAnew);
	cout << "size(JA): " << result[1].size()  << endl;
	displayVector(JAnew);
	result.push_back(V);
	cout << "size(V): " << result[2].size()  << endl;
	displayVector(V);
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
		if(i%1000 == 0) cout << "y[" << i << "]= " << y[i] << endl;
	}
}

int main(int argc, char *argv[])
{
	string filename = "test.mtx"; // "s3rmt3m1.mtx"
	vector<vector<double>> vec = readMatrix(filename);
	int vec_size = vec[1].size() - 1;
	vector<double> x(vec_size, 1.0);
	vector<double> y(vec_size, 0.0);
	if (filename == "test.mtx"){
		matrixmulCCS(vec[2],vec[0],vec[1],x,y);
	}
	else
	{
		v_matrixmulCCS(vec[2],vec[0],vec[1],x,y);
	}

    return 0;
}