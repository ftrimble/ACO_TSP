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
KRATOS=-O7 -DKRATOS -lm
DEBUG=$(KRATOS) -DDEBUG_MODE
BLUE=-O3 -DBLUE

pdfs:
	pdflatex findings; bibtex findings; pdflatex findings; pdflatex findings

all:
	mpicc $(SOURCES) -o $(OUTPUT) $(KRATOS)
# run with "mpirun -np 4 ./aco_tsp input/FILE_NAME NUM_CITIES"

debug: 
	mpicc $(SOURCES) $(DEBUG) -o $(OUTPUT)

run:
	mpirun -np $(RANKS) ./aco_tsp input/$(INPUT).tsp $(NCITIES)

path_dist:
	gcc src/tsp_path_distance.c -o tsp_path_distance -Wall -lm

run_path_dist:
	./tsp_path_distance input/$(INPUT).tsp input/$(INPUT).opt.tour $(NCITIES)

blue:
	mpicc $(SOURCES) $(BLUE) -o $(BLUEOUTDIR)/$(OUTPUT)

clean:
	rm -r *~ $(OUTPUT)
