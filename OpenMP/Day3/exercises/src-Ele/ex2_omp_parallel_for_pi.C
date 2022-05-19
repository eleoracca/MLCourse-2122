//---------------------------------------------------------------------
// This program is made to compute pi numerical value.
// It will numerically compute the integral of
//
//                   4/(1+x*x)
//
// from 0 to 1 using variable number of steps.
// The value of this integral is pi. The knowing of pi is quite 
// precise and the result can be compared to the "real" value.
//
// This parallel version uses omp parallel for.
// It uses a timer class defined in a separate header file.
//
// History: 
// Written by Gabriele Gaetano Fronzé, 01/18.
// Based on Tim Mattson implementation, 11/99.
//---------------------------------------------------------------------

#include "StopWatch.h"
#include <omp.h>
#include <iostream>
#include "common.h"

int main()
{
    StopWatch stopWatch;

    double pi, sum, s = 0.0;
    double step = 1.0/(double) common::num_steps; //x-step
    int n_threads=1;


    #pragma omp parallel
    {
        n_threads = omp_get_num_threads();
        uint thread_id = omp_get_thread_num();
    /*
        // OpenMP can handle for us for loop ranges!
        // TIP: you have to add something in the following line...

        #pragma omp for
        for (unsigned long long i=0; i<=common::num_steps; i++) {
            double x = (i - 0.5) * step; //computing the x value
            sum += 4.0 / (1.0 + x * x); //adding to the cumulus
        }
    */

    /*
        #pragma omp for reduction (+:sum)
        for (unsigned long long i=0; i<=common::num_steps; i++) {            
            double x = (i - 0.5) * step; //computing the x value
            sum += 4.0 / (1.0 + x * x); //adding to the cumulus
        }    
    */
        
    /*
    To avoid raise conditions, we can use also update and critical, but it takes a lot more time, as you need to wait for all the threads to access the variable one per time.
    Reduction though, creates a sort of vector of the sums and gives one item to all each of the threads. At the end it will sum all 
    
    more steps -> enlarges the time non linearly
    it is because of the sequential part of the code, which is still there
    */
        
    /*
        #pragma omp for
        for (unsigned long long i=0; i<=common::num_steps; i++) {
            double x = (i - 0.5) * step; //computing the x value
            s = 4.0 / (1.0 + x * x);
            
            #pragma omp atomic update
            sum += s; //adding to the cumulus
        }
    */
        
    ///*
        #pragma omp for
        for (unsigned long long i=0; i<=common::num_steps; i++) {
            double x = (i - 0.5) * step; //computing the x value
            s = 4.0 / (1.0 + x * x);
            
            #pragma omp critical
            sum += s; //adding to the cumulus
        }
    //*/
    
    }

    pi = step * sum;

    common::print_results(pi, n_threads);
    
    return 0;
    
}
