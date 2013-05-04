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

#define ITER_MAX 10000000
#define IMPROVE_REQ 1000000

#ifdef KRATOS
double clock_rate = 2666700000.0; 
#else /* Using Blue Gene/Q */
double clock_rate = 1600000000.0; 
#endif

struct city;

double **distances = NULL;
double **inverted_distances = NULL;
double **pheromones = NULL;
double **pheromones_recv = NULL; // secondary matrix for holding the allreduce

int num_cities;
struct city * cities;

double rho, alpha, beta;
  
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

// after calling this function, distances and inverted_distances will be filled
// returns -1 on failure and 0 on success
int parse_input_file(char* input_file_path) {
  int i, j;
  char * this_line;
  char buffer[100];
    
  FILE *input_file = fopen(input_file_path, "rb");
  if (input_file == NULL) {
    perror ("The following error occurred");
    return -1;
  }
    
  // get header information
  this_line = fgets(buffer, 100, input_file); // NAME
  this_line = fgets(buffer, 100, input_file); // TYPE
  this_line = fgets(buffer, 100, input_file); // COMMENT
  this_line = fgets(buffer, 100, input_file); // DIMENSION
  sscanf(this_line,"DIMENSION: %d",&num_cities);
  this_line = fgets(buffer, 100, input_file); // EDGE_WEIGHT_TYPE
  this_line = fgets(buffer, 100, input_file); // NODE_COORD_SECTION

  cities = (struct city*)calloc(num_cities, sizeof(struct city));

  // read cities
  for (i=0; i < num_cities; i++) {
    this_line = fgets(buffer, 100, input_file);
        
    sscanf(this_line, " %d %lf %lf", 
	   &(cities[i].id_number), 
	   &(cities[i].x_coord), 
	   &(cities[i].y_coord)
	   );
  }

  fclose(input_file);
  return 0;
}


/* |                                                         | */
/* ----------------------------------------------------------- */


// this finds the probability that an ant will travel along the edge
// between loc and dest, given the locations it has already visited.
double edge_prob(int loc, int dest, int *visited) { 
  int i; 
  double num, denom = 0; 
  
  num = pow(pheromones[loc][dest],alpha)*inverted_distances[loc][dest]; 
  for ( i = 0; i < num_cities; ++i ) { 
    if ( visited[i] || i == loc ) continue; 
    denom += pow(pheromones[loc][i],alpha)*inverted_distances[loc][i];
  } 

  return num/denom; 
}

// finds the next edge to visit
int findEdge(int loc, int *visited) {
  double prob = genrand_res53(),    // rolls the dice
    totProb = 0;    
  int y;                            // next city to go to

  for ( y = 0; y < num_cities; ++y ) {
    if ( y == loc || visited[y] ) continue; // does not solve TSP

    // follows an edge if it corresponds with our dice roll.
    double toProb = edge_prob(loc,y, visited); 
    if ( toProb > prob ) return y;
    prob -= toProb;
    totProb += toProb;
  }

  // this should never be reached because the probabilities are 
  // normalized to sum to 1. Nonetheless, if such a thing does
  // occur, a notification will be useful.
  return -1; 
}

double findPath(int *visited, int *path) {
  int curr_loc,                      // current location
    city_start,                      // starting location
    i = 0;
  int num_to_visit;                  // cities left
  double dist = 0;                   // path length

  curr_loc = city_start = 
    genrand_int32() % num_cities;    // start ant randomly
  num_to_visit = num_cities -1;      // don't include current city

  path[i++] = city_start;
  
  // reset visited; only passed to avoid repeated memory allocation
  memset(visited,0,num_cities*sizeof(int));
  visited[city_start] = 1;

  // finds a TSP path
  while ( num_to_visit ) {
    int dest = findEdge(curr_loc,visited); // next location
    if ( dest == -1 ) {
      fprintf(stderr,"ERROR: Could not find a suitable edge\n");
      return -1;
    }

    // updates path length and pheromones
    dist += distances[curr_loc][dest];
    pheromones[dest][curr_loc] = 
      pheromones[curr_loc][dest] += inverted_distances[curr_loc][dest];
        
    // moves to next city and checks it off
    visited[dest] = 1;
    curr_loc = dest;
    --num_to_visit;

    path[i++] = curr_loc;
  }
  
  // loops back to starting city
  dist += distances[curr_loc][city_start];
  pheromones[curr_loc][city_start] = 
    pheromones[city_start][curr_loc] += inverted_distances[curr_loc][city_start];
  
  return dist;
}


// this function finds the distance of the "optimal" TSP
double findTSP(int taskid ) {
  int num_iters = 0,                                       // number of attempts
    i, j,                                                  // counters
    time_since_improve = 0,                                // time since last improvement
    *visited = (int *)calloc(num_cities,sizeof(int)),      // visited cities
    *path = (int *)calloc(num_cities,sizeof(int)),
    *optPath = (int *)calloc(num_cities,sizeof(int)),
    sendOptPath;
  double dist,                                             // path length
    best_distance = DBL_MAX,                               // best path
    last_improve = best_distance;                          // last path improvement

  while ( num_iters < ITER_MAX ) {
    // pheromones decay 
    for ( i = 0; i < num_cities; ++i )
      for ( j = 0; j < num_cities; ++j )
	pheromones[i][j] = rho*pheromones[i][j];
      
    dist = findPath(visited, path);
    if ( dist == -1 ) dist = best_distance;

    // updates min distance over all processors
    if ( dist >= best_distance ) {
      sendOptPath = 0;
      dist = best_distance;
    }
    else  sendOptPath = 1;
    MPI_Allreduce(&dist,&best_distance,1,MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // checks for end conditions
    if ( best_distance < last_improve ) {
      last_improve = best_distance;
      time_since_improve = 0;
    }
    else ++time_since_improve;
    if ( time_since_improve > IMPROVE_REQ ) return best_distance;

    // grabs the new pheromones over all processors
    MPI_Allreduce(&pheromones[0][0],&pheromones_recv[0][0],num_cities*num_cities,
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // received pheromones are our new pheromone values
    double **tmp = pheromones_recv;
    pheromones_recv = pheromones;
    pheromones = tmp;
      
    ++num_iters;
  }

  return best_distance;
}

int main(int argc, char *argv[]) {

  /* --------------------------------------------------------------------- */
  /* | Variable declaration and initialization                           | */
    
  MPI_Status status;          // status for MPI_Test
  int taskid, numtasks;       // MPI info
  int input_parse_return_val; // error check on input parsing
  int i, j;                   // loop indeces
  double best_distance;       // used for termination and output
  double execTime;            // used to time program execution

  // Random Number Seeding
  unsigned long rng_init_seeds[6] = {0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
  unsigned long rng_init_length = 6;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  execTime = rdtsc();

  if (argc != 2) {
    if ( taskid == 0 )
      fprintf(stderr, "ERROR: Invalid usage.\nUSAGE: %s INPUT\n", argv[0]);

    return EXIT_FAILURE;
  }

  alpha = beta = 1;
  rho  = 1.0/numtasks;

  // seed MT19937
  rng_init_seeds[0] = taskid;
  init_by_array(rng_init_seeds, rng_init_length);
    
  /* |                                                                   | */
  /* --------------------------------------------------------------------- */
    
  /* --------------------------------------------------------------------- */
  /* | Read input file, initialize distances, and broadcast              | */
    
  if (taskid == 0) input_parse_return_val = parse_input_file(argv[1]);
  MPI_Barrier(MPI_COMM_WORLD);

  // Check for errors (couldn't open the input file)
  MPI_Bcast(&input_parse_return_val, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (input_parse_return_val != 0) { MPI_Finalize(); exit(EXIT_FAILURE); }

  // gets cities info across all processors
  MPI_Bcast(&num_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if ( taskid != 0 ) cities = (struct city *)calloc(num_cities,sizeof(struct city));
  MPI_Bcast(cities, num_cities*sizeof(struct city), MPI_BYTE, 0, MPI_COMM_WORLD);

  // allocates memory for our matrices
  distances = alloc2dcontiguous(num_cities, num_cities);
  inverted_distances = alloc2dcontiguous(num_cities, num_cities);
  pheromones = alloc2dcontiguous(num_cities, num_cities);
  pheromones_recv = alloc2dcontiguous(num_cities, num_cities);
    
  // calculate and store distances, their inverses and pheromones.
  for (i=0; i < num_cities; i++) {
    for (j=0; j < num_cities; j++) {
      distances[i][j] = city_distance(cities[i], cities[j]);
      inverted_distances[i][j] = i == j ? 0 : pow(1/distances[i][j],beta);
      pheromones[i][j] = pheromones_recv[i][j] = 1;
    }
  }

  free(cities);

  // prints out the matrix that was broadcasted in DEBUG_MODE
#ifdef DEBUG_MODE 
  if (taskid == 0) {
    printf("distances:\n");
    print_matrix(distances, num_cities, num_cities);
    printf("\n");
  }
#endif
    
  /* |                                                                   | */
  /* --------------------------------------------------------------------- */
    
  /* --------------------------------------------------------------------- */
  /* | Algorithm processing                                              | */

  best_distance = findTSP(taskid);
  if (taskid == 0) printf("distance: %f\n", best_distance);
    
  /* |                                                                   | */
  /* --------------------------------------------------------------------- */

  /* --------------------------------------------------------------------- */
  /* | Cleanup and output                                                | */
    
  // calculates the times in seconds that we used.
  execTime = (rdtsc() - execTime)/clock_rate;
    
  // rank 0 performs output
  if (taskid == 0)
    printf("Finished, program took %f seconds.\n", execTime);
    
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
