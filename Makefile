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
REPORTOUT=findings.pdf
REPORT=findings.aux \
       findings.bbl \
       findings.blg \
       findings.log \
       findings.pdf

all: $(OUTPUT)
quick: kratos
kratos: $(OUTPUT) $(REPORTOUT)

aco_tsp:
	mpicc $(SOURCES) -o $(OUTPUT) $(KRATOS)


findings.pdf:
	./run

debug: 
	mpicc $(SOURCES) $(DEBUG) -o $(OUTPUT)

path_dist:
	gcc src/tsp_path_distance.c -o tsp_path_distance -Wall -lm

run_path_dist:
	./tsp_path_distance input/$(INPUT).tsp input/$(INPUT).opt.tour $(NCITIES)

blue:
	mpicc $(SOURCES) $(BLUE) -o $(BLUEOUTDIR)/$(OUTPUT)

clean:
	rm -r $(REPORT) *~ $(OUTPUT)

