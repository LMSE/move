# MoVE

MoVE is a tool that identifies metabolic valves to switch from a growth to a production stage. It makes use of constrained minimal cut sets.
Associated research article: Venayak, N., von Kamp, A., Klamt, S. et al. MoVE identifies metabolic valves to switch between phenotypic states. Nat Commun 9, 5332 (2018). https://doi.org/10.1038/s41467-018-07719-4

# Setup

CellNetAnalyzer requires additions to both java and matlab path to run. Firstly, in the script CNAPATH\startcna.m (CNAPATH is the directory where your CNA code is) the two following variables need to be set accordingly:

cplex_matlab_path= 'CPLEXROOT\cplex\matlab\OS'; % CPLEXROOT is the directory where you have installed cplex, OS varies depending on your computer version
cplex_jar_file= 'CPLEXROOT\cplex\lib\cplex.jar';

This startcna.m script has to be run whenever matlab is started. It is also run as part of the MoVE_startup.m script.

Secondly, you have to add the cplex shared library to java.library.path. First locate the file librarypath.txt shown below.

	MATLAB\version\toolbox\local\librarypath.txt	(MATLAB is the directory where Matlab is installed)
	
Then make a copy of this text file to your desktop and open it. Add a line at the end of the file:

	CPLEXROOT\cplex\bin\OS	

Then copy this file back to its original location and override the original. This change will not take place until Matlab is restarted.

To test if the connections are properly setup, run setup_cplex_inner_class_access() (after having called starcna). If no errors are thrown, then your paths have been properly setup.

Finally, set the paths in MoVE_startup.m according to your installation and run this script to start up everything.


# update - March 2020
geneVa extends the existing approach for the computation of metabolic valves. It uses gene-product-reaction associations to compute static knockouts and valves at the gene-level, instead of the reaction-level. (https://github.com/LMSE/geneVa) - Author: Philipp Schneider - https://github.com/VonAlphaBisZulu 
