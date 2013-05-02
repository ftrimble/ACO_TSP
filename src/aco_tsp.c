/******************************************************************************/
/******************************************************************************/
/***** Forest Trimble,               ******************************************/
/***** Scott Todd,                   ******************************************/
/***** {trimbf,todds}@rpi.edu        ******************************************/
/***** Project:                      ******************************************/
/*****   Ant Colony Optimization,    ******************************************/
/*****   Travelling Salesman Problem ******************************************/
/***** Due: May 7, 2013              ******************************************/
/******************************************************************************/
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // requires -lm at compilation
#include <unistd.h>
#include <float.h>
#include <mpi.h>

#include "rdtsc.h"
#include "MT19937.h"

#ifdef KRATOS
double clock_rate = 2666700000.0; 
#else /* Using Blue Gene/Q */
double clock_rate = 1600000000.0; 
#endif

double **distances = NULL;
double **inverted_distances = NULL;
double **pheromones = NULL;
double **pheromones_recv = NULL;

int num_cities;

double rho = 0.5, 
  alpha = 1, beta = 1;


/* ----------------------------------------------------------- */
/* | Distance setup tools                                    | */
struct city {
    int id_number;
    double x_coord, y_coord;
};

double city_distance(struct city city1, struct city city2) {
    return sqrt( pow((city1.x_coord - city2.x_coord), 2) + 
                 pow((city1.y_coord - city2.y_coord), 2) );
}

double invert_double(double num) {
    if (num > 0.0001) { return 1 / num; }
    else { return num; }
}
/* |                                                         | */
/* ----------------------------------------------------------- */


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

void print_matrix(double ** matrix, int rows, int cols) {
    int i, j;
    for (i=0; i < rows; i++) {
        for (j=0; j < cols; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

double edge_prob(int loc, int dest, int *visited) { 
   // uses alpha = beta = 1 
   int i; 
   double num, denom; 
   num = pheromones[loc][dest]*inverted_distances[loc][dest]; 
   for ( i = 0; i < num_cities; ++i ) { 
     if ( visited[i] || i == loc ) continue; // not interested in already visited locations 
     denom += pheromones[loc][i]*inverted_distances[loc][i];
   } 

   return num/denom; 
 }

// finds the next edge to visit
int findEdge(int prob, int loc, int *visited) {
  int y;
  for ( y = 0; y < num_cities; ++y ) {
    if ( y == loc || visited[y] ) continue;
    double toProb = edge_prob(loc,y, visited);
    if ( toProb > prob ) return y;
    prob -= toProb;
  }
}


// after calling this function, distances and inverted_distances will be filled
// returns -1 on failure and 0 on success
int parse_input_file(char* input_file_path, int num_cities) {
    int i, j;
    char * this_line;
    char buffer[100];
    
    struct city* cities = (struct city*)calloc(num_cities, sizeof(struct city));
    
    FILE *input_file = fopen(input_file_path, "rb");
    if (input_file == NULL) {
        perror ("The following error occurred");
        return -1;
    }
    
    // skip header information
    this_line = fgets(buffer, 100, input_file); // NAME
    this_line = fgets(buffer, 100, input_file); // TYPE
    this_line = fgets(buffer, 100, input_file); // COMMENT
    this_line = fgets(buffer, 100, input_file); // DIMENSION
    this_line = fgets(buffer, 100, input_file); // EDGE_WEIGHT_TYPE
    this_line = fgets(buffer, 100, input_file); // NODE_COORD_SECTION
    
    // read cities
    for (i=0; i < num_cities; i++) {
        this_line = fgets(buffer, 100, input_file);
        
        sscanf(this_line, "%d %lf %lf", 
              &(cities[i].id_number), 
              &(cities[i].x_coord), 
              &(cities[i].y_coord)
             );
    }
    
    // calculate and store distances and their inverses
    for (i=0; i < num_cities; i++) {
        for (j=0; j < num_cities; j++) {
            distances[i][j] = city_distance(cities[i], cities[j]);
            inverted_distances[i][j] = invert_double(distances[i][j]);
        }
    }
    
    free(cities);
    fclose(input_file);
    
    return 0;
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        fprintf(stderr, "ERROR: Invalid usage.\n");
        fprintf(stderr, "USAGE: %s input_file num_cities\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* --------------------------------------------------------------------- */
    /* | Variable declaration and initialization                           | */
    
    MPI_Status status; // status for MPI_Test
    int taskid, numtasks; // MPI info
    int input_parse_return_val;
    int i, j; // loop indeces
        
    // int num_cities; made global
    sscanf(argv[2], "%d", &num_cities);
    
    double best_distance = DBL_MAX; // used for termination and output
    
    double execTime; // used to time program execution
    
    // random number seeding
    unsigned long rng_init_seeds[6] = {0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
    unsigned long rng_init_length = 6;
    
    distances = alloc2dcontiguous(num_cities, num_cities);
    inverted_distances = alloc2dcontiguous(num_cities, num_cities);
    pheromones = alloc2dcontiguous(num_cities, num_cities);
    pheromones_recv = alloc2dcontiguous(num_cities, num_cities);
    
    // initialize pheromone graph
    // distances and inverted distances are initialized by rank 0 and broadcast
    for (i=0; i < num_cities; i++) {
        for (j=0; j < num_cities; j++) {
            pheromones[i][j] = 1;
        }
    }
    
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
    init_by_array(rng_init_seeds, rng_init_length);
    
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------- */
    /* | Read input file, initialize distances, and broadcast              | */
    
    if (taskid == 0) {
        input_parse_return_val = parse_input_file(argv[1], num_cities);
    }
    
    // Check for errors (couldn't open the input file)
    MPI_Bcast(&input_parse_return_val, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (input_parse_return_val != 0) { MPI_Finalize(); exit(EXIT_FAILURE); }
    
    MPI_Bcast(&(distances[0][0]), num_cities * num_cities, MPI_DOUBLE, 
              0, MPI_COMM_WORLD);
    MPI_Bcast(&(inverted_distances[0][0]), num_cities * num_cities, MPI_DOUBLE, 
              0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // DEBUG
    // if (taskid != 0) { print_matrix(distances, num_cities, num_cities); }
    // printf("task #%d\n", taskid);
    // print_matrix(distances, num_cities, num_cities);
    
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */
    
    /* --------------------------------------------------------------------- */
    /* | Algorithm processing                                              | */
    
    
    int num_iters = 0, 
      ITER_MAX = 1000; // small test value
    while ( num_iters < ITER_MAX ) {
      int city_start = genrand_int32() % num_cities, num_to_visit = num_cities; 
      double dist = 0;
      int *visited = (int *)calloc(sizeof(int),num_cities);
      int i, j;
      for ( i = 0; i < num_cities; ++i )
        for ( j = 0; j < num_cities; ++j )
          pheromones[i][j] = (1 - rho)*pheromones[i][j];
      
      
      // finds a path
      while ( num_to_visit ) {
        double rand = genrand_res53(); // rolls the dice
        
        // next location
        int dest = findEdge(rand, city_start,visited);

        // updates path length and pheromones
        dist += distances[city_start][dest];
        pheromones[city_start][dest] += inverted_distances[city_start][dest];
        
        // moves to next city and checks it off
        visited[dest] = 1;
        city_start = dest;
        --num_to_visit;
      }

      dist = dist < best_distance ? dist : best_distance;

      // grabs the new pheromones over all processors
      MPI_Allreduce(&pheromones[0][0],&pheromones_recv[0][0],num_cities*num_cities,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      double **tmp = pheromones_recv;
      pheromones_recv = pheromones;
      pheromones = tmp;

      // updates min distance over all processors
      MPI_Allreduce(&dist,&best_distance,1,MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      
      ++num_iters;
    }
    
    
    
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */

    /* --------------------------------------------------------------------- */
    /* | Cleanup and output                                                | */
    
    // calculates the times in seconds that we used.
    execTime = (rdtsc() - execTime)/clock_rate;
    
    // rank 0 performs output
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
     
    free(&distances[0][0]);                     free(distances);
    free(&inverted_distances[0][0]);            free(inverted_distances);
    free(&pheromones[0][0]);                    free(pheromones);
    
    
    MPI_Finalize();

    return EXIT_SUCCESS;
    /* |                                                                   | */
    /* --------------------------------------------------------------------- */
}
