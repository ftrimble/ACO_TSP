set terminal post eps solid

set xlabel 'Number of Tasks'
set logscale x 2

set title 'Distance vs. Number of Tasks in a Strong Scaling Study'
set ylabel 'Distance'
set output 'plots/data_dist.eps'
plot 'plots/data_dist.dat' w l title ''

set logscale y 2
set title 'Time vs. Number of Tasks in a Strong Scaling Study'
set ylabel 'Time (Seconds)'
set output 'plots/data_time.eps'
plot 'plots/data_time.dat' lc rgb 'red'   w l title    'Overall Average',    \
     'plots/data_comm.dat' lc rgb 'blue'  w l title 'Communication Average', \
     'plots/data_comp.dat' lc rgb 'green' w l title  'Computation Average'  

##########################################################
##### Plots some unbordered graphs to show the tours #####
##########################################################

unset logscale xy
set border 0
set xlabel ''
set ylabel ''
set title ''
unset ytics
unset xtics

set output 'plots/pr2392.eps'
plot 'plots/pr2392.dat' w l title ''

set title 'The optimal tour for a 150 city problem'
set output 'plots/opt_tour.eps'
plot 'plots/opt_tour.dat' w l title ''

set title 'The ACO tour for a 150 city problem'
set output 'plots/data_tour.eps'
plot 'plots/data_tour.dat' w l title ''

##########################################################
##########################################################
##########################################################