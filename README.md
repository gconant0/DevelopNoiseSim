# DevelopNoiseSim
Simulations of organismal development in the presence of expression noise.

This program simulates the development of a multicellular organism of  2^*b(t+1)-1) cells, given the user specifying t (the number of cell-type transitions) and b (the number of cell divisions between cell-type transitions). A prefect development is first simulated. Then a noisy gene expression model is used to allow failure of cell transitions. The Stochkit package is used for these simulations (see below). A Postscript visualization of the development is also made. 

Command line usage:
embryo_sim_wrand <T> <B> <ERRORMODEL.XML> <ERRORRATE>

where 0< ERRORRATE <1

Alternatively, 

embryo_sim_wrand <T> <B> <ERRORMODEL.XML> -n:<DF#>

where DF# is the cutoff for the diff. factor that causes differentiation failure.

DEPENDANCIES:
Stochkit: https://github.com/StochSS/StochKit
(Assumes that tau_leaping_exp_adapt_serial is in your path)

GNU PlotUtils: https://www.gnu.org/software/plotutils/

Compiling:
Example compile command with g++:
g++ -o embryo_sim_wrand embryo_sim_wrand.cpp -lplot
