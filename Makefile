###############################################################################
###############################################################################
##### Forest Trimble,                   #######################################
##### Scott Todd                        #######################################
##### {trimbf,todds}@rpi.edu            #######################################
##### Project:                          #######################################
#####   Ant Colony Optimization,        #######################################
#####   Travelling Salesman Problem     #######################################
##### Due: May 7, 2013                  #######################################
###############################################################################
###############################################################################

SOURCES=src/aco_tsp.c
OUTPUT=aco_tsp
BLUEOUTDIR=~/data-sb
KRATOS=-O7 -DKRATOS -lm -Wall -pthread
DEBUG=$(KRATOS) -DDEBUG_MODE
BLUE=-O3 -DBLUE

all:
	mpicc $(SOURCES) -o $(OUTPUT) $(KRATOS)
# run with "mpirun -np 4 ./aco_tsp input/FILE_NAME NUM_CITIES"

debug: 
	mpicc $(SOURCES) $(DEBUG) -o $(OUTPUT)

run_test:
	mpirun -np $(RANKS) ./aco_tsp input/$(INPUT).tsp

path_dist:
	gcc src/tsp_path_distance.c -o tsp_path_distance -Wall -lm

run_path_dist:
	./tsp_path_distance input/$(INPUT).tsp input/$(INPUT).opt.tour $(NCITIES)

blue:
	mpicc $(SOURCES) $(BLUE) -o $(BLUEOUTDIR)/$(OUTPUT)

clean:
	rm -r        \
	findings.bbl \
	findings.aux \
	findings.blg \
	findings.log \
	findings.pdf \
	*~ $(OUTPUT)

