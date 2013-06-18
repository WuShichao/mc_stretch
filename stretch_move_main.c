/*
Copyright (c) 2013, Alex Kaiser
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer. Redistributions
 in binary form must reproduce the above copyright notice, this list
 of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "stretch_move_sampler.h"
#include "cl-helper.h"
#include "stretch_move_util.h"
#include "constants.h"
#include "initial_conds.h"

//    Stretch Move MCMC sampler in OpenCL
//    Alex Kaiser, Courant Institute, 2012
//    Email: user: adkaiser  domain: gmail

 

int main(int argc, char **argv){


    // User set parameters
    cl_int chain_length = 200000;                     // Allocate to store this much chain, sampler runs this many steps at once

    int burn_length = 1000000;                        // Length of burn in

    cl_int annealing_loops = 20;                      // Run this many temperatures of simulated annealing
    cl_int steps_per_loop = 50000;                    // This many steps each

    cl_int dimension = N_TH + N_STEPS * NX ;          // Dimension of the state vector
    cl_int walkers_per_group = 1024;                  // Total number of walkers is twice this
    size_t work_group_size = 32;                      // Work group size. Use 1 for CPU, larger number for GPU
    double a = 1.4;                                   // Coefficient for range of 'z' random variable
    cl_int pdf_number = 0;                            // Does not matter in this example
    const char *plat_name = CHOOSE_INTERACTIVELY;
    const char *dev_name  = CHOOSE_INTERACTIVELY;


    // set this parameter for a debug run
    int easy = 0;
    if(easy){
        chain_length     = 10000;
        burn_length      = 10000;
        steps_per_loop   = 1000;
    }



    // set parameters about which components to save
    cl_int num_to_save        = N_TH;
    cl_int *indices_to_save   = (cl_int *) malloc(num_to_save * sizeof(cl_int));
    for(int i=0; i<num_to_save; i++)
        indices_to_save[i] = i;


    // read the observations from a file
    // note: N_OBS, NY are defined in "constants.h"
    cl_int data_length = N_OBS * NY;
    cl_float *data = (cl_float *) malloc(data_length * sizeof(cl_float)) ;
    if(!data){ perror("Allocation failure obs_temp"); abort(); }
    char file_name[] = "noisy_data.txt";
    read_arrays(data, N_OBS, NY, file_name);


    // Initialize the sampler
    sampler *samp = initialize_sampler(chain_length, dimension, walkers_per_group, work_group_size, a, pdf_number,
                                       data_length, data, num_to_save, indices_to_save, plat_name, dev_name);


    // initialize the initial conditions in the struct
    // initial conditions are written in constand array in definitions
    for(int i=0; i<NX; i++)
        (samp->data_st)->x_initial[i] = X_INITIAL_DEF[i] ;


    // Run simulated annealing too speed up convergence
    // set a generic cooling schedule, {1/10, 1/9 ... 1}
    cl_float *cooling_schedule = (cl_float *) malloc(annealing_loops * sizeof(cl_float));
    int idx=0;
    for(int i=annealing_loops; i>0; i--) cooling_schedule[idx++] = 1.0f / ( (cl_float) i);

    // run the annealing
    run_simulated_annealing(samp, cooling_schedule, annealing_loops, steps_per_loop);
    free(cooling_schedule);



    run_burn_in(samp, burn_length);


    // run the sampler
    run_sampler(samp);

    // --------------------------------------------------------------------------
    // The array samp->samples_host now contains samples ready for use.
    //
    // Array is in component major order.
    // To access the i-th saved sample of sample j use
    //     samp->samples_host[i + j*(samp->num_to_save)]
    //
    // Dimension is (samp->N x samp->total_samples)
    // --------------------------------------------------------------------------


    // print summary of the run including basic some statistics
    print_run_summary(samp);

    // run acor to estimate autocorrelation time
    run_acor(samp);

    // Output some histograms to Matlab, don't output gnuplot
    char matlab_hist = 1, gnuplot_hist = 0;
    output_histograms(samp, matlab_hist, gnuplot_hist);

    // free resources
    free_sampler(samp);

    return 1;
}



