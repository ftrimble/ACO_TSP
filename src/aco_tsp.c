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
#include <pthread.h>
#include <mpi.h>

#include "rdtsc.h"
#include "MT19937.h"

int ITER_MAX = 256;
int IMPROVE_REQ = 256;
#define num_threads 8
#define SUCCESS 1
#define FAIL 0

#ifdef KRATOS
double clock_rate = 2666700000.0; 
#else /* Using Blue Gene/Q */
double clock_rate = 1600000000.0; 
#endif

// mostly, this argument is in place to avoid repeated
// memory allocation and freeing. the dist argument
// actually serves as a return value, while the remaining
// arguments are arrays that would otherwise need to be
// allocated and freed at each iteration. 
struct TSPargs {
  double dist;
  double **nums;
  double *denoms;
  int *visited;
  int *path;
};

struct CompArg {
  int row,
    doRows;
  double rho;
};

struct city {
  int id_number;
  double x_coord, y_coord;
};

unsigned long long startTime;

// data
int num_cities;
double **distances = NULL;
double **inverted_distances = NULL;
double **nums;
double **pheromones = NULL;
double **pheromones_recv = NULL; // secondary matrix for holding the allreduce
pthread_mutex_t **pher_mutex;
struct city * cities;

// heuristics
double rho, rhothreads, rhoprocs, 
  alpha, 
  beta;

// mpi stuff
int taskid, numtasks;

// best solutions
double best_distance, my_best_distance;
int *my_best_path;

/* ----------------------------------------------------------- */
/* | Miscellaneous Utility Functions                         | */


/* This function allocates a contiguous block of memory *
 * so that MPI_Isend/Irecv operations are much easier.  */
double ** alloc2dcontiguous(int rows, int cols) {
  int i;
  double *data = (double *)calloc(rows*cols, sizeof(double));
  double **array = (double **)calloc(rows, sizeof(double*));

  for (i = 0; i < rows; ++i)
    array[i] = &data[cols*i];

  return array;
}

void print_matrix(double ** matrix, int rows, int cols) {
  int i, j;
  printf("********************** MATRIX ******************************\n");
  for (i=0; i < rows; i++) {
    for (j=0; j < cols; j++) {
      printf("%e ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n***********************************************************\n\n");
}

void print_path(int * path) {
  int i;
  for(i=0; i < num_cities; ++i) {
    // if results don't line up, check for off-by-one here
    printf("%f %f\n", cities[path[i]].x_coord, cities[path[i]].y_coord);
  }
  printf("%f %f\n", cities[path[0]].x_coord, cities[path[0]].y_coord);
    
}

/* |                                                         | */
/* ----------------------------------------------------------- */


/* ----------------------------------------------------------- */
/* | Distance setup tools                                    | */

double city_distance(struct city city1, struct city city2) {
  return sqrt( pow((city1.x_coord - city2.x_coord), 2) + 
	       pow((city1.y_coord - city2.y_coord), 2) );
}

// after calling this function, distances and inverted_distances will be filled
// returns -1 on failure and 0 on success
int parse_input_file(char* input_file_path) {
  int i;
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

/* ----------------------------------------------------------- */
/* | TSP/ACO Functions                                       | */


// this function symmetrically adds pherQ 
// to the correct location in the pheromone matrix.
// simplifies using mutex locks
void addPheromones(int x, int y, double pherQ) {
  pthread_mutex_lock(&pher_mutex[x][y]);
  pheromones[x][y] += pherQ;
  pthread_mutex_unlock(&pher_mutex[x][y]);

  pthread_mutex_lock(&pher_mutex[y][x]);
  pheromones[y][x] += pherQ;
  pthread_mutex_unlock(&pher_mutex[y][x]);
}

// decays the pheromones
void *decayPheromones(void *args) {
  int i, j;
  struct CompArg *arg = (struct CompArg *)args;
  
  for ( i = arg->row; i < arg->row + arg->doRows; ++i )
    for ( j = 0; j < num_cities; ++j )
      pheromones[i][j] = arg->rho*pheromones[i][j];

return NULL;
}

// gets numerators of the probabilities.
void *getFracts(void *args) {
  int i, j;
  struct CompArg *arg = (struct CompArg *)args;

  /*************************************************
   *             tau_ij^alpha*d(C_i,C_j)^beta      *
   * p_ij = -------------------------------------- *
   *        sum_n=1^m tau_in^alpha*d(C_i,C_n)^beta *
   *************************************************/
  for ( i = arg->row; i < arg->row + arg->doRows; ++i )
    for ( j = 0; j < num_cities; ++j )
      nums[i][j] = pow(pheromones[i][j],alpha)*pow(inverted_distances[i][j],beta);

  return NULL;
}

// finds the next edge to visit
int findEdge(int loc, double denom, struct TSPargs *arg) {
  double prob = genrand_res53(),    // rolls the dice
    totProb = 0;    
  int i;
  
  for ( i = 0; i < num_cities; ++i ) {
    if ( i == loc || arg->visited[i] ) continue; // does not solve TSP
    
    // follows an edge if it corresponds with our dice roll.
    double toProb = nums[loc][i]/denom;

    if ( isnan(toProb) ) {
      return -2;
      break;
    } 

    if ( toProb > prob )
      return i;

    prob -= toProb;
    totProb += toProb;
  }
  
  // this should never be reached because the probabilities are 
  // normalized to sum to 1. Nonetheless, if such a thing does
  // occur, a notification will be useful.
  return -1;
}

void * findPath(void *args) {
  int curr_loc,                          // current location
    city_start,                          // starting location
    i = 0, j;
  int num_to_visit;                      // cities left
  struct TSPargs *arg = (struct TSPargs *)args;
  double denom;

  // resets the values of the argument
  memset(arg->visited,0,num_cities*sizeof(int));
  arg->dist = 0;

  curr_loc = city_start = 
    genrand_int32() % num_cities;    // start ant randomly
  num_to_visit = num_cities -1;      // don't include current city

  arg->visited[city_start] = 1;
  arg->path[i++] = city_start;

  // finds a TSP path
  while ( num_to_visit ) {
    denom = 0;
    for ( j = 0; j < num_cities; ++j ) 
      if ( !arg->visited[j] ) denom += nums[curr_loc][j];
    
    int dest = findEdge(curr_loc, denom, arg); // next location

    if ( dest == -1 ) {
      fprintf(stderr,"ERROR: Could not find a suitable edge\n");
      arg->dist = -1;
      return NULL;
    }
    if ( dest == -2 ) {
      fprintf(stderr,"%s %s", 
	      "Pheromones shrinking too rapidly,", 
	      "rho needs adjusting. Exiting...\n");
      arg->dist = -2;
      return NULL;
    }
    
    // updates path length and pheromones
    arg->dist += distances[curr_loc][dest];
    addPheromones(curr_loc, dest, inverted_distances[curr_loc][dest]);
        
    // moves to next city and checks it off
    arg->path[i++] = dest;
    arg->visited[dest] = 1;
    curr_loc = dest;
    --num_to_visit;
  }
  
  // loops back to starting city
  arg->dist += distances[curr_loc][city_start];
  addPheromones(curr_loc, city_start, inverted_distances[curr_loc][city_start]);

  return NULL;
}


// this function finds the distance of the best path found before termination
void findTSP(int* num_iters, unsigned long long* total_communication_cycles, 
               unsigned long long* total_computation_cycles) {

  int i, j,                                  // counters
    time_since_improve = 0;                  // time since last improvement
  double last_improve = best_distance;       // last path improvement
  unsigned long long startComputeCycles,
    startCommCycles,
    endCommCycles;
  pthread_t *threads;
  struct TSPargs *path_args;
  struct CompArg *compute_args;

  // memory for threads
  threads = (pthread_t *)calloc(num_threads, sizeof(pthread_t));
  path_args = (struct TSPargs *)calloc(num_threads, sizeof(struct TSPargs));
  compute_args = (struct CompArg *)calloc(num_threads,sizeof(struct CompArg));

  for ( i = 0; i < num_threads; ++i ) {
    path_args[i].path = (int *) calloc(num_cities,sizeof(int));
    path_args[i].visited = (int *)calloc(num_cities,sizeof(int));
    compute_args[i].doRows = num_cities/num_threads;
    compute_args[i].row = num_cities/num_threads*i;
    if ( i == num_threads - 1 ) compute_args[i].doRows += num_cities % num_threads;
  }

  while ( *num_iters < ITER_MAX ) {

    // ------- Timing ------- //
    startComputeCycles = rdtsc();

    for ( i = 0; i < num_threads; ++i )
      compute_args[i].rho = rhothreads;

    // individual threaded sub-colony problem attempt
    for ( i = 0; i < IMPROVE_REQ; ++i ) {

      // pheromones decay 
      for ( j = 0; j < num_threads; ++j ) 
	pthread_create(&threads[j], NULL, &decayPheromones, 
		       (void *)&compute_args[j]);
      for ( j = 0; j < num_threads; ++j ) 
	pthread_join(threads[j], NULL);

      // calculate nums and denoms
      for ( j = 0; j < num_threads; ++j ) 
	pthread_create(&threads[j], NULL, &getFracts, 
		       (void *)&compute_args[j]);
      for ( j = 0; j < num_threads; ++j )
	pthread_join(threads[j], NULL);

      // each thread finds a path 
      for ( j = 0; j < num_threads; ++j ) 
	pthread_create(&threads[j], NULL, &findPath,
		       (void *)&path_args[j]);

      for ( j = 0; j < num_threads; ++j ) {
	pthread_join(threads[j],NULL);
	if ( path_args[j].dist == -2 ) {
	  *num_iters = ITER_MAX;
	  break;
	}
	else if ( path_args[j].dist < my_best_distance && 
		  path_args[j].dist != -1 ) {
	  my_best_distance = path_args[j].dist;
	  memcpy(my_best_path,path_args[j].path,num_cities*sizeof(int));
	}
      }
    }

    // ------- Timing ------- //
    *total_computation_cycles += 
      rdtsc() - startComputeCycles;
    startCommCycles = rdtsc();
    // ------- Timing ------- //

#ifdef DEBUG_MODE
    /* --------------- Debugging output ----------------- */
    if (taskid == 0)
      printf("%s%d %s %f %s%f\n",
	     "Before MPI_Allreduce: Iteration #",*num_iters,
	     "time:",(startCommCycles-startTime)/clock_rate,
	     "best_distance:",best_distance);
    /* --------------- Debugging output ----------------- */
#endif // DEBUG_MODE

    // updates min distance over all processors
    MPI_Allreduce(&my_best_distance,&best_distance,1,
		  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // checks for end conditions
    if ( best_distance < last_improve ) {
      last_improve = best_distance;
      time_since_improve = 0;
    }
    else ++time_since_improve;
    if ( time_since_improve > IMPROVE_REQ ) return;

    // grabs the new pheromones over all processors
    MPI_Allreduce(&pheromones[0][0],&pheromones_recv[0][0],num_cities*num_cities,
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // received pheromones are our new pheromone values
    double **tmp = pheromones_recv;
    pheromones_recv = pheromones;
    pheromones = tmp;

    // ------- Timing ------- //
    endCommCycles = rdtsc();
    *total_communication_cycles += 
      endCommCycles - startCommCycles;
    // ------- Timing ------- //

    // pheromones decay 
    for ( i = 0; i < num_threads; ++i ) {
      compute_args[i].rho = rhoprocs;
      pthread_create(&threads[i], NULL, 
		     &decayPheromones, 
		     &compute_args[i]);
    }

#ifdef DEBUG_MODE      
    /* --------------- Debugging output ----------------- */
    if (taskid == 0) {
      printf("%s%d %s %f %s %f\n",
	     "After MPI_Allreduce: Iteration #",*num_iters,
	     "time:",(endCommCycles-startTime)/clock_rate,
	     "best_distance:",best_distance);
        
      for(i=0; i < num_cities; ++i) printf("%d ", my_best_path[i]);
      print_matrix(pheromones, num_cities, num_cities);
      printf("\n");
    }
    /* --------------- Debugging output ----------------- */
#endif // DEBUG_MODE
    
    ++(*num_iters);
  }
}

/* |                                                         | */
/* ----------------------------------------------------------- */


int main(int argc, char *argv[]) {

  /* --------------------------------------------------------------------- */
  /* | Variable declaration and initialization                           | */
    
  MPI_Status status;          // status for MPI_Test
  int taskid, numtasks,       // MPI info
    input_parse_return_val,   // error check on input parsing
    i, j,                     // loop indeces
    signal;
  double execTime;            // used to time program execution

  best_distance = my_best_distance = DBL_MAX;
  
  // Random Number Seeding
  unsigned long rng_init_seeds[6] = {0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
  unsigned long rng_init_length = 6;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

  startTime = execTime = rdtsc();

  if (argc != 2) {
    if ( taskid == 0 )
      fprintf(stderr, "ERROR: Invalid usage.\nUSAGE: %s INPUT\n", argv[0]);

    return EXIT_FAILURE;
  }
  
  // strong scaling study
  ITER_MAX /= numtasks;
  IMPROVE_REQ = ITER_MAX;

  // heuristic constants
  alpha = 1;
  beta = 16;

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
  nums = alloc2dcontiguous(num_cities, num_cities);

  pher_mutex = (pthread_mutex_t **)calloc(num_cities,sizeof(pthread_mutex_t*));
  for ( i = 0; i < num_cities; ++i )
    pher_mutex[i] = (pthread_mutex_t *)calloc(num_cities,sizeof(pthread_mutex_t));

  double max_id, min_id = DBL_MAX, avg_id = 0, id;
  // calculate and store distances, their inverses and pheromones.
  for (i=0; i < num_cities; i++) {
    for (j=0; j < num_cities; j++) {
      distances[i][j] = city_distance(cities[i], cities[j]);
      
      id = i == j ? 0 : 1/distances[i][j];
      avg_id += inverted_distances[i][j] = id;
      if ( id < min_id ) min_id = id;
      if ( id > max_id ) max_id = id;
      
      pthread_mutex_init(&pher_mutex[i][j], NULL);
    }
  }
  
  for ( i = 0; i < num_cities; ++i ) 
    for ( j = 0; j < num_cities; ++j ) 
    pheromones[i][j] = pheromones_recv[i][j] = avg_id;

  
  rho = min_id;
  rhothreads = avg_id/(avg_id+num_threads*rho); 
  rhoprocs = avg_id/(avg_id+numtasks*rho);

  /* |                                                                   | */
  /* --------------------------------------------------------------------- */
    
  /* --------------------------------------------------------------------- */
  /* | Algorithm processing                                              | */

  my_best_path = (int *)calloc(num_cities,sizeof(int));
  int num_iters = 0;
  
  unsigned long long total_communication_cycles = 0; // MPI_Allreduce time
  unsigned long long total_computation_cycles = 0;   // computation time

  findTSP(&num_iters, &total_communication_cycles, 
	  &total_computation_cycles);

  /* |                                                                   | */
  /* --------------------------------------------------------------------- */

  /* --------------------------------------------------------------------- */
  /* | Cleanup and output                                                | */
  
  double my_communication_time = total_communication_cycles/clock_rate;
  double my_computation_time = total_computation_cycles/clock_rate;
  
  double avg_communication_time; // across all MPI_Ranks
  double avg_computation_time;   // across all MPI_Ranks
  
  MPI_Reduce( &my_communication_time,
	      &avg_communication_time,
	      1,
	      MPI_DOUBLE,
	      MPI_SUM,
	      0,
	      MPI_COMM_WORLD
	      );
  MPI_Reduce( &my_computation_time,
	      &avg_computation_time,
	      1,
	      MPI_DOUBLE,
	      MPI_SUM,
	      0,
	      MPI_COMM_WORLD
	      );
  
  if (taskid == 0) {
    execTime = (rdtsc() - execTime)/clock_rate;
    printf("%d\n", num_iters);
    printf("%f\n", execTime);
    printf("%f\n", avg_communication_time/numtasks);
    printf("%f\n", avg_computation_time/numtasks);
    printf("%f\n", best_distance);
  }
  
  // determine which rank has the best path and output that path
  
  // each non-rank 0 worker check with rank 0 if it has the optimal tour
  // this is to avoid duplicate optimums being output multiple times
  if(taskid != 0) {
    // send my best distance to rank 0
    MPI_Send(&my_best_distance, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
    // wait for a message from rank 0, 0 == I'm not the best, 1 == I'm the best
    MPI_Recv(&signal, 1, MPI_INT, 0, 
	     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
             
    if ( signal ) print_path(my_best_path);
  }
  
  // rank 0 performs output and coordinates who has the optimal tour
  if (taskid == 0) {
    int already_sent = 0;
    double recv_distance;
    signal = FAIL;
    
    if (my_best_distance == best_distance ) {
      print_path(my_best_path);
      already_sent = 1;
    }
    
    for(i=1; i < numtasks; ++i) {
      MPI_Recv(&recv_distance, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD,
	       &status);
        
      if (recv_distance == best_distance && !already_sent)
	signal = already_sent = SUCCESS;
      else signal = FAIL;
      MPI_Send(&signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  
  /* |                                                                   | */
  /* --------------------------------------------------------------------- */

    
  /* --------------------------------------------------------------------- */
  /* | Memory deallocation and exit                                      | */
    
  free(cities);
  free(my_best_path);
  for ( i = 0; i < num_cities; ++i ) {
    for ( j = 0; j < num_cities; ++j ) 
      pthread_mutex_destroy(&pher_mutex[i][j]);
    free(pher_mutex[i]);
  }
  free(pher_mutex);
  
  /* frees up the memory allocated for our arrays. *
   * recall that the array was initiated in one    *
   * contiguous chunk, so one call to free should  *
   * deallocate the whole underlying structure.    */
  free(&nums[0][0]);                      free(nums);
  free(&distances[0][0]);                 free(distances);
  free(&inverted_distances[0][0]);        free(inverted_distances);
  free(&pheromones[0][0]);                free(pheromones);
  free(&pheromones_recv[0][0]);           free(pheromones_recv);
    
  MPI_Finalize();

  return EXIT_SUCCESS;

  /* |                                                                   | */
  /* --------------------------------------------------------------------- */
}
