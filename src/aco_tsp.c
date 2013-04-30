/*******************************************************************************/
/*******************************************************************************/
/***** Forest Trimble,               *******************************************/
/***** Scott Todd,                   *******************************************/
/***** {trimbf,todds}@rpi.edu        *******************************************/
/***** Project:                      *******************************************/
/*****   Ant Colony Optimization,    *******************************************/
/*****   Travelling Salesman Problem *******************************************/
/***** Due: May 7, 2013              *******************************************/
/*******************************************************************************/
/*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <pthread.h>

#include "rdtsc.h"
#include "MT19937.h"

#ifdef KRATOS
double clock_rate = 2666700000.0; 
#else /* Using Blue Gene/Q */
double clock_rate = 1600000000.0; 
#endif

double **distances = NULL;
double **inverted_distances = NULL;


/* This function allocates a contiguous block of memory *
 * so that MPI_Isend/Irecv operations are much easier.  */
double ** alloc2dcontiguous(int rows, int cols) {
    int i;
    double *data = (double *)calloc(rows*cols, sizeof(double));
    double **array = (double **)calloc(rows, sizeof(double));

    for (i = 0; i < rows; ++i) {
        array[i] = &data[cols*i];
    }

    return array;
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
        fprintf(stderr, "ERROR: Invalid usage.\n");
        fprintf(stderr, "USAGE: %s input_file\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* --------------------------------------------------------------------- */
    /* | Variable declaration and initialization                           | */
    
    int taskid, numtasks; // MPI info
    int i, j; // loop indeces
    
    double execTime; // used to time program execution
    
    // random number seeding
    unsigned long rng_init_seeds[6] = {0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
    unsigned long rng_init_length = 6;
    

    /* |                                                                   | */
    /* --------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------- */
    /* | MPI and RNG initialization                                        | */
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    execTime = rdtsc();

    // seed MT19937
    rng_init_seeds[0] = taskid;
    init_by_array(rng_init_seeds,rng_init_length);
    
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------- */
    /* | Read input file, initialize matrices, and broadcast               | */
    
    if (taskid == 0) {
        char input_file_path[100];
        strcpy(input_file_path, argv[1]);
        
    }
    
    
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------- */
    /* | Algorithm processing                                              | */
    
    
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */

    /* --------------------------------------------------------------------- */
    /* | Cleanup and output                                                | */
    
    
    // calculates the times in seconds that we used.
    execTime = (rdtsc() - execTime)/clock_rate;
    
    // node 0 performs output
    if (taskid == 0) {
        printf("Finished, program took %f seconds.\n", execTime);
    }
    
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */

    
    /* --------------------------------------------------------------------- */
    /* | Memory deallocation and exit                                      | */
    
    /* frees up the memory allocated for our arrays. *
     * recall that the array was initiated in one    *
     * contiguous chunk, so one call to free should  *
     * deallocate the whole underlying structure.    */
     
    // free(&distances[0][0]);                     free(distances);
    // free(&inverted_distances[0][0]);            free(inverted_distances);
    free(distances);
    free(inverted_distances);
    
    
    MPI_Finalize();

    return EXIT_SUCCESS;
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */
}
