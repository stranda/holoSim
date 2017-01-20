# holoSim
R scripts for coupling spatially explicit forward simulations with coalescent simulations

#Contents
holoSim-Ash.R
	This is the main one.  Start here.  It sources in other scripts and runs nrep simulations.  
	This script can be run from the command line with, for instance...
		'Rscript holoSim-Ash.R 1 100'
	in this case, the output would be Ash_1.out and the script would attempt to run 100 replicates
	the first argument is the "node"... change this to avoid overwriting results
	the second argument is the number of replicates to attempt
	The script writes files line by line, so 'wc -l Ash_1.out' will give you an idea of how many successful runs you've done (header row gets counted too, so subtract one for the # of replicates)
	This script also saves an R data object for each replicate (Ash_<node>-<replicate>.Rda) that contains the gtypes object for the simulation (genotypes, strataG format) and the parameters used in the replicate
	‘Rsource’ and ‘SIMS’ objects point to folders where the scripts are and where you want the simulation files to be written, respectively.  
	#Change those before running anything!
	
colrate.R
	10 different options for colonization rate, all are calculated in holoSim-Ash.R

empty_or_sampled.R
	has vectors for 
		empty cells (Atlantic ocean, outside of Fraxinus range, etc.)
		grid cell IDs for sampled populations
		sample sizes for each sampled population
		
holoStats.R
	function to calculate stats from the simulation output -- strataG object
	should work with observed Fraxinus data too
		segregating sites per pop
		private segregating sites per pop & total
		observed heterozygosity per pop & total
		pairwise Fst
		total Fst
	Also contains a function to mask simulated data to match the pattern of missing data for Fraxinus.  This depends on a ‘data mask’ .csv file that has TRUE/FALSE where the observed data is missing/present, respectively.
			
runFSC.R
	function that uses strataG functions to write a fastsimcoal file, based on the pops object from rmetasim, includes option to specify a growth model — exponential growth or step change in population size

	
getpophist3.R
	minor modifications from Allan's getpophist script.  Differences here are...
		includes colonization rate window, measures abundance at beginning and end of this period
		new logical test for when to stop simulating (empty cells will have is.na(pops$arrive) == TRUE at the end of the sim)
	
make-landscape2.R
	minor modifications from Allan's makelandscape script.  Differences here are...
		increased fecundity to keep populations around, 1.2 now, previously 1.05
		landscape.popcoords() fxn edited to have popIDs increase from left to right and bottom to top
			bottom row, lower left corner -> lower right corner -> second row left side -> etc.
	
plothist.R
	the original from Allan's work.  Not changed, included so you can check out colonization history, if interested.
	Note, it doesn't like NA's in the pops object
		Use: plothist(pops[is.na(pops$arrive) == FALSE,]) 
