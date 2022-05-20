#include <omp.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <vector>
#include <numeric>

using namespace std;
using namespace std::chrono;


// -----------------------------------------------------------------------------
// -------------------- Declarations
// -----------------------------------------------------------------------------

double factorial(int n); // computes the factorial

void range(int max, int n_max); // parallel executions with ranges for each thread
void forsimple(int max, int n_max); // parallel execitions with pragma for

void forcritical(int max, int n_max); // parallel executions with pragma for and critical statement
void foratomic(int max, int n_max); // parallel executions with pragma for and atomic statement
void forreduction(int max, int n_max); // parallel executions with pragma for and reduction statement


// -----------------------------------------------------------------------------
// -------------------- Implementations
// -----------------------------------------------------------------------------

int main(){

	// maximum values
	int max = 1; // number of threads
	int n_max = 100000; // iterations of the sum
	
	#pragma omp parallel
	{
		#pragma omp single
		max = omp_get_num_threads();
	}
	
	// comparison between for loop divided by the user and the pragma for
	range(max, n_max);
	forsimple(max, n_max);
	
	// assessing the race conditions
	forcritical(max, n_max);
	foratomic(max, n_max);
	forreduction(max, n_max);
	
	cout << "\n -----------------------------------------------------------------" << endl;
	cout << "          Execution completed" << endl;
	return 24;
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// -------------------- Factorial function
double factorial(int n){
	if(n == 1){
		return 1;
	}
	else{
		return n * factorial(n-1);
	}
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// -------------------- Range function
void range(int max, int n_max){
	cout << "\n -----------------------------------------------------------------" << endl;
	cout << "          Starting range function" << endl;
	
	// output file for the timing 
	FILE *fp;
	char file[32];
	sprintf (file, "./Output/range.dat");
	fp = fopen(file,"w");
	
	// variables for timing
	high_resolution_clock::time_point begin;
	high_resolution_clock::time_point end;
	double times[max];
	
	// support variables
	vector<double> e;
	e.resize(1);
	double euler;
	
	uint thread_id = omp_get_thread_num();
	
	// for loop on the number of threads
	for(int j = 1; j <= max; j++){
	
		// resetting the variables
		fill(e.begin(), e.end(), 0.0); // reset the sum at each cycle
		begin = high_resolution_clock::now(); // the time begins
		
		#pragma omp parallel num_threads(j)
		{	
			#pragma omp single
			e.resize(j);
			
			thread_id = omp_get_thread_num();
			
			for(int i = thread_id; i < n_max; i += j){
				e[thread_id] += 1.0/factorial(i+1);
			}
		
		}
		
		euler = accumulate(e.begin(), e.end(), 1.0); // adding the sums to get Napier's value	
		end = high_resolution_clock::now(); // the time ends
		
		// computing elapsed time
		duration<double, std::milli> temp = end - begin;
		times[j] = temp.count() / 1000.;
	
		fprintf(fp, "%12.6e %d \n", times[j], j);
	
		cout << "Elapsed time: " << setprecision(7) << times[j] << " s" << endl;
		cout << "Euler's constant value = " << euler << endl << endl;
	}
	
	fclose(fp);

}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// -------------------- Forsimple function
void forsimple(int max, int n_max){
	cout << "\n -----------------------------------------------------------------" << endl;
	cout << "          Starting forsimple function" << endl;
	
	// output file for the timing 
	FILE *fp;
	char file[32];
	sprintf (file, "./Output/forsimple.dat");
	fp = fopen(file,"w");
	
	// variables for timing
	high_resolution_clock::time_point begin;
	high_resolution_clock::time_point end;
	double times[max];
	
	// support variables
	vector<double> e;
	double euler;
	
	uint thread_id = omp_get_thread_num();
	
	// for loop on the number of threads
	for(int j = 1; j <= max; j++){
	
		// resetting the variables
		fill(e.begin(), e.end(), 0.0); // reset the sum at each cycle
		begin = high_resolution_clock::now(); // the time begins
		
		#pragma omp parallel num_threads(j)
		{	
			#pragma omp single
			e.resize(j);
			
			thread_id = omp_get_thread_num();
			
			#pragma omp for
			for(int i = 1; i <= n_max; i++){
				e[thread_id] += 1.0/factorial(i);
			}
			
		}
		
		euler = accumulate(e.begin(), e.end(), 1.0); // adding the sums to get Napier's value	
		end = high_resolution_clock::now(); // the time ends
		
		// computing elapsed time
		duration<double, std::milli> temp = end - begin;
		times[j] = temp.count() / 1000.;
	
		fprintf(fp, "%12.6e %d \n", times[j], j);
	
		cout << "Elapsed time: " << setprecision(7) << times[j] << " s" << endl;
		cout << "Euler's constant value = " << euler << endl << endl;
	}
	
	fclose(fp);
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// -------------------- Forcritical function
void forcritical(int max, int n_max){
	cout << "\n -----------------------------------------------------------------" << endl;
	cout << "          Starting forcritical function" << endl;
	
	// output file for the timing 
	FILE *fp;
	char file[32];
	sprintf (file, "./Output/forcritical.dat");
	fp = fopen(file,"w");
	
	// variables for timing
	high_resolution_clock::time_point begin;
	high_resolution_clock::time_point end;
	double times[max];
	
	// support variables
	double e = 0.0;
	double es = 0.0;
	double euler;
	
	// for loop on the number of threads
	for(int j = 1; j <= max; j++){
	
		// resetting the variables
		e = 0.0; // reset the sum at each cycle
		es = 0.0; // reset the sum at each cycle
		begin = high_resolution_clock::now(); // the time begins
		
		#pragma omp parallel num_threads(j)
		{	
			#pragma omp for
			for(int i = 1; i <= n_max; i++){
				e = 1.0/factorial(i);
				
				#pragma omp critical
				es += e;
			}
		
		}
		
		euler = es + 1.0; // adding the first addendum of the series, that is 1/0! (= 1)		
		end = high_resolution_clock::now(); // the time ends
		
		// computing elapsed time
		duration<double, std::milli> temp = end - begin;
		times[j] = temp.count() / 1000.;
	
		fprintf(fp, "%12.6e %d \n", times[j], j);
	
		cout << "Elapsed time: " << setprecision(7) << times[j] << " s" << endl;
		cout << "Euler's constant value = " << euler << endl << endl;
	}
	
	fclose(fp);
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// -------------------- Foratomic function
void foratomic(int max, int n_max){
	cout << "\n -----------------------------------------------------------------" << endl;
	cout << "          Starting foratomic function" << endl;
	
	// output file for the timing 
	FILE *fp;
	char file[32];
	sprintf (file, "./Output/foratomic.dat");
	fp = fopen(file,"w");
	
	// variables for timing
	high_resolution_clock::time_point begin;
	high_resolution_clock::time_point end;
	double times[max];
	
	// support variables
	double e = 0.0;
	double es = 0.0;
	double euler;
	
	// for loop on the number of threads
	for(int j = 1; j <= max; j++){
	
		// resetting the variables
		e = 0.0; // reset the sum at each cycle
		es = 0.0; // reset the sum at each cycle
		begin = high_resolution_clock::now(); // the time begins
		
		#pragma omp parallel num_threads(j)
		{	
			#pragma omp for
			for(int i = 1; i <= n_max; i++){
				e = 1.0/factorial(i);
				
				#pragma omp atomic update
				es += e;
			}
		
		}
		
		euler = es + 1.0; // adding the first addendum of the series, that is 1/0! (= 1)		
		end = high_resolution_clock::now(); // the time ends
		
		// computing elapsed time
		duration<double, std::milli> temp = end - begin;
		times[j] = temp.count() / 1000.;
	
		fprintf(fp, "%12.6e %d \n", times[j], j);
	
		cout << "Elapsed time: " << setprecision(7) << times[j] << " s" << endl;
		cout << "Euler's constant value = " << euler << endl << endl;
	}
	
	fclose(fp);
}
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// -------------------- Forreduction function
void forreduction(int max, int n_max){
	cout << "\n -----------------------------------------------------------------" << endl;
	cout << "          Starting forreductin function" << endl;
	
	// output file for the timing 
	FILE *fp;
	char file[32];
	sprintf (file, "./Output/forreduction.dat");
	fp = fopen(file,"w");
	
	// variables for timing
	high_resolution_clock::time_point begin;
	high_resolution_clock::time_point end;
	double times[max];
	
	// support variables
	double e = 0.0;
	double euler;
	
	// for loop on the number of threads
	for(int j = 1; j <= max; j++){
	
		// resetting the variables
		e = 0.0; // reset the sum at each cycle
		begin = high_resolution_clock::now(); // the time begins
		
		#pragma omp parallel num_threads(j)
		{	
			#pragma omp for reduction(+:e) 
			for(int i = 1; i <= n_max; i++){
				e += 1.0/factorial(i);
			}
				
		}
		
		euler = e + 1.0; // adding the first addendum of the series, that is 1/0! (= 1)
		end = high_resolution_clock::now(); // the time ends
		
		// computing elapsed time
		duration<double, std::milli> temp = end - begin;
		times[j] = temp.count() / 1000.;
	
		fprintf(fp, "%12.6e %d \n", times[j], j);
	
		cout << "Elapsed time: " << setprecision(7) << times[j] << " s" << endl;
		cout << "Euler's constant value = " << euler << endl << endl;
	}
	
	fclose(fp);
}

