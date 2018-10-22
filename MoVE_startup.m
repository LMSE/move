%%% Modify the paths accordingly for your installation
%%% CNA: https://www2.mpi-magdeburg.mpg.de/projects/cna/download.html
%%% Ensure paths are also set correctly in CellNetAnalyzer/startcna.m lines 5 and 6
cna_root = ''; % e.g. 'CellNetAnalyzer'
cobra_toolbox_path = ''; %e.g. cobra_2.0.5\cobra
MoVE_path = fullfile(pwd(), 'MoVE');

addpath(genpath(MoVE_path));

cd(cobra_toolbox_path);
initCobraToolbox();

cd(cna_root);
startcna(1);

cd(MoVE_path);
changeCobraSolver('glpk');
