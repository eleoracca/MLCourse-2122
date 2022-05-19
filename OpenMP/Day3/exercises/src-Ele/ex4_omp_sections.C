#include "StopWatch.h"
#include <omp.h>
#include <iostream>
#include <thread>
#include <chrono>

using namespace std::chrono;

#define MAX omp_get_max_threads()

int main(){
    StopWatch stopwatch;
    
    const long maxIteration = 100000;
    bool debug_output=false;

    // This is needed to allow the presence of a parallel for nested in another parallel region
    omp_set_nested(1);

    int outputCounter = 0;
    bool done = false;

    double pi, sum, s = 0.0;
    double step = 1.0/(double) maxIteration; //x-step
    int n_threads=1;

    #pragma omp parallel sections shared(outputCounter,done) num_threads(3)
    {
        
        #pragma omp section
        {
            n_threads = omp_get_num_threads();
            uint thread_id = omp_get_thread_num();
            
            if ( omp_get_thread_num()==0 ){ 
                if(debug_output) printf("monitor th_id: %d, n_th: %d\n",omp_get_thread_num(), omp_get_num_threads());
                while(!done){
                    printf("%d/%ld\r",outputCounter,maxIteration);
                    std::this_thread::sleep_for (microseconds(50));
                }
                printf("%d/%ld\n",outputCounter,maxIteration);
            }
            
        }

        #pragma omp section
        {
            #pragma omp parallel num_threads(MAX)
            {
                if(debug_output) printf("computing th_id: %d, n_th: %d\n",omp_get_thread_num(), omp_get_num_threads());
                #pragma omp for
                for (unsigned long long i=0; i<=maxIteration; i++) {
                    double x = (i - 0.5) * step; //computing the x value
                    s = 4.0 / (1.0 + x * x);
            
                    #pragma omp atomic update
                    sum += s; //adding to the cumulus

                    #pragma omp atomic update
                    outputCounter++;
                }
            }
            done = true;
        }
    }

    return 0;
}