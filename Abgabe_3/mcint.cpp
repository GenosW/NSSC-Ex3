#include <assert.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <omp.h>

using namespace std;


int main(int argc, char *argv[]) 
{
	struct timespec ts, te;
	double runtime;
	int threads = 1;
    double sum = 0;
    double c = 0;
	vector<double> T(threads);

	// parse command line arguments
	assert(argc == 7);
	string impl = "UNKNOWN";
	double Xmax = -1;
	double Xmin = -1;
	double Ymax = -1;
	double Ymin = -1;
	double sampltemp = -1;
	int tid;
	omp_set_num_threads(threads);
	{
		istringstream tmp(argv[1]);
		tmp >> impl;
	}
	{
		istringstream tmp(argv[2]);
		tmp >> Xmin;
    }
	{
		istringstream tmp(argv[3]);
		tmp >> Xmax;
	}
	{
		istringstream tmp(argv[4]);
		tmp >> Ymin;
	}
	{
		istringstream tmp(argv[5]);
		tmp >> Ymax;
	}
	{
		istringstream tmp(argv[6]);
		tmp >> sampltemp;
	}
	int sampl = (int) sampltemp; 
	cout << "impl=" << impl << endl;
	cout << "Xmin=" << Xmin << endl;
	cout << "Xmax=" << Xmax << endl;
	cout << "Ymin=" << Ymin << endl;
	cout << "Ymax=" << Ymax << endl;
	cout << "sampl=" << sampl << endl;
	assert(sampl > 0);
	

	if (impl == "SINX")
	{
		mt19937_64 MASRNG(1000);
		for(int i=0; i<threads; ++i) 
		{
			 T[i] = MASRNG();
		}
		c = (Xmax-Xmin)/sampl;
		clock_gettime(CLOCK_REALTIME, &ts);
		    #pragma omp parallel private(tid) shared(sampl,Xmax,Xmin,Ymax,Ymin) reduction(+:sum)
			{		
				tid = omp_get_thread_num();
				mt19937_64 THRRNG(T.at(tid));
				uniform_real_distribution<double> dice(0.0, 1.0);
					for(int i=0; i<sampl; ++i) 
					{
						double zw = 0;
						double temps = sin((Xmax-Xmin)*dice(THRRNG)+Xmin);
						if(temps > Ymin && temps < Ymax)
						{
							zw = temps;
						}
						if(temps < Ymin)
						{
							zw = Ymin;
						}
						if(temps > Ymax)
						{
							zw = Ymax;
						}
						sum = sum + zw;
					}
			cout << "Erg by threads SINX: " << c*sum << endl;				
			}
			cout << "Total Erg SINX: " << c*sum/threads << endl;
			clock_gettime(CLOCK_REALTIME, &te);		
	}
	else if (impl == "SIN2XINV")
	{
		mt19937_64 MASRNG(1000);
		for(int i=0; i<threads; ++i) 
		{
			 T[i] = MASRNG();
		}
		c = (Xmax-Xmin)/sampl;
		clock_gettime(CLOCK_REALTIME, &ts);
		    #pragma omp parallel private(tid) shared(sampl,Xmax,Xmin,Ymax,Ymin) reduction(+:sum)
			{		
				tid = omp_get_thread_num();
				mt19937_64 THRRNG(T.at(tid));
				uniform_real_distribution<double> dice(0.0, 1.0);
					for(int i=0; i<sampl; ++i) 
					{
						double zw = 0;
						double temps = sin(1/((Xmax-Xmin)*dice(THRRNG)+Xmin))*sin(1/((Xmax-Xmin)*dice(THRRNG)+Xmin));
						if(temps > Ymin && temps < Ymax)
						{
							zw = temps;
						}
						if(temps < Ymin)
						{
							zw = Ymin;
						}
						if(temps > Ymax)
						{
							zw = Ymax;
						}
						sum = sum + zw;
					}
			cout << "Erg by threads SIN2XINV: " << c*sum << endl;				
			}
			cout << "Total Erg SIN2XINV: " << c*sum/threads << endl;
			clock_gettime(CLOCK_REALTIME, &te);		
	}
	else if (impl == "XCUBE")
	{
		mt19937_64 MASRNG(1000);
		for(int i=0; i<threads; ++i) 
		{
			 T[i] = MASRNG();
		}

		c = (Xmax-Xmin)/sampl;
		clock_gettime(CLOCK_REALTIME, &ts);
		    #pragma omp parallel private(tid) shared(sampl,Xmax,Xmin,Ymax,Ymin) reduction(+:sum)
			{		
				tid = omp_get_thread_num();
				mt19937_64 THRRNG(T.at(tid));
				uniform_real_distribution<double> dice(0.0, 1.0);
					for(int i=0; i<sampl; ++i) 
					{
						double zw = 0;
						double temps = ((Xmax-Xmin)*dice(THRRNG)+Xmin)*((Xmax-Xmin)*dice(THRRNG)+Xmin)*((Xmax-Xmin)*dice(THRRNG)+Xmin);
						if(temps > Ymin && temps < Ymax)
						if(temps > Ymin && temps < Ymax)
						{
							zw = temps;
						}
						if(temps < Ymin)
						{
							zw = Ymin;
						}
						if(temps > Ymax)
						{
							zw = Ymax;
						}
						sum = sum + zw;
					}
			cout << "Erg by threads XCUBE: " << c*sum << endl;				
			}
			cout << "Total Erg XCUBE: " << c*sum/threads << endl;
			clock_gettime(CLOCK_REALTIME, &te);
	}
	else
	{
		cout << "Implementation not recognized. Terminating!";
		return -1;
	}
	runtime = double(te.tv_sec - ts.tv_sec) + double(te.tv_nsec - ts.tv_nsec)/1e9;
	cout << "Implementation: " << impl << endl;
	cout << "Runtime of Integration: " << runtime << endl;
	cout << "Number of Threads: " << threads << endl;
	return 0;
}


