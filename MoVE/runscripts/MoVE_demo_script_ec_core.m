% % Core model
load ecoli_core
model = changeRxnBounds(model,'EX_glc(e)',-10,'l');
reac_off = [];
del_exchanges = [];

target = 'EX_succ(e)'
timeout = Inf;	% s

solution_type = 3; % Search for feasible solutions. Optimal solution is returned if infinite timeout provided

% This defines the type of algo to run, can explicitely define a number of
% valves, or max number of valves (unweighted), or choose a number of
% valves based on a weight. Should probably use weighted first to see if
% any solution exists
objective_type = 'unweighted';  % 'weighted' or 'unweighted'

    % if weighted
    objective_weight = 1;

    % if unweighted
    number_of_valves = 3;
    valve_constraint = 'less_than_or_equal'; % 'equal' or 'less_than_or_equal'

growth_stage_min_mue = 0.5;	% percent of max growth rate

% use a low production_stage_min_yield when checking 
% whether coupling is possible at all
production_stage_min_yield = 0.5; % percent of max growth rate
production_stage_min_mue = 0; % biomass yield (gdw/mol)

% Number of threads to use, 0 for auto
num_threads = 0;

% Plot the production envelope of the result
% Note no visible envelope will be displayed if mu = 0
plot_envelope = true;

% Array of reaction indices which should not be used as valves. Empty string for none
no_valve = '';

% Alteratively a single_valve id can be chosen, only this valve will be allowed
single_valve = '';	% e.g. 'EX_o2(e)'

% Extra reactions not to knockout, can still be valves
no_KO = []; % e.g. {'EX_o2(e)'}

% Calculate and return solution (this is solved via JVM)
total_solution = calc_MoVE(model, target, 'EX_glc(e)', growth_stage_min_mue, ...
    production_stage_min_yield, production_stage_min_mue, objective_type, ...
    objective_weight, number_of_valves, valve_constraint, 3, ...
    timeout, num_threads, single_valve, no_valve, plot_envelope, false,...
    reac_off, del_exchanges, no_KO);