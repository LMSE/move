function [ total_solution ] = calc_MoVE( model, target, glucose, growth_stage_min_mue, ...
    production_stage_min_yield, production_stage_min_mue, objective_type, ...
    objective_weight, number_of_valves, valve_constraint, solution_type, timeout,...
    num_threads, single_valve_reaction, no_valve2, plot_envelope, distributed, ...
    reac_off, del_exchanges, no_KO)
    % MoVE Calculates multi-objective cutsets for a single product
    % Error codes
    % 0: Preparation error
    % 1: FVA error
    % 2: LP error
    % 3: Parameter value error
    % 4: No solution found from MILP
    % 9: Other error

    %% This is the configuration information for cplex's distributed search.
    %% These variables must be set for a distributed search
    if distributed
        node_1 = getenv('NODE_1')
        node_2 = getenv('NODE_2')
        node_3 = getenv('NODE_3')
        node_4 = getenv('NODE_4')

        vmconfig = sprintf(['<?xml version="1.0" encoding="US-ASCII"?>'...
                    '<vmc>'...
                    '  <machine name="machine1">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '  <machine name="machine2">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '  <machine name="machine3">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '  <machine name="machine4">'...
                    '    <transport type="TCP/IP">'...
                    '      <address host="%s" port="8050" />'...
                    '    </transport>'...
                    '  </machine>'...
                    '</vmc>'],node_1,node_2,node_3,node_4);
    end

    total_solution = [];
    total_solution.solution = 0;
    total_solution.name = target;

    %% Set cplex parameters
    cplex_inner= setup_cplex_inner_class_access();
    emphasizeAccuracy= cplex_inner.ParameterSet.constructor.newInstance([]);
    emphasizeAccuracy.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
    % emphasizeAccuracy.setParam(cplex_inner.BooleanParam.PreInd, false);
    emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpOpt, 1e-9);
    emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpRHS, 1e-9);
    % emphasizeAccuracy.setParam(cplex_inner.IntParam.ScaInd, 1);
    emphasizeAccuracy.setParam(cplex_inner.IntParam.AdvInd, 0);

    % explicit parameter set for solving the MILP
    MILPparameters= cplex_inner.ParameterSet.constructor.newInstance([]);
    MILPparameters.setParam(cplex_inner.IntParam.MIPEmphasis, cplex_inner.MIPEmphasis.Balanced);
    MILPparameters.setParam(cplex_inner.IntParam.MIPDisplay, 2);
    MILPparameters.setParam(cplex_inner.DoubleParam.EpInt, 0);
    MILPparameters.setParam(cplex_inner.DoubleParam.EpRHS, 1e-9);
    MILPparameters.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
    MILPparameters.setParam(cplex_inner.DoubleParam.WorkMem, 1024*10);
    MILPparameters.setParam(cplex_inner.IntParam.RandomSeed, 2);
    MILPparameters.setParam(cplex_inner.IntParam.Threads,num_threads); % limit threads for testing
    MILPparameters.setParam(cplex_inner.IntParam.NodeFileInd, 3); % Saves nodefile to disk (2 for uncompressed, 3 for compressed)
    if distributed
        MILPparameters.setParam(cplex_inner.DistMIP.Rampup.TimeLimit, timeout/4);
    end
    num_evs_lim = 5; %max number of solutions

    %% Prepare the model
    biomassRxnName = model.rxns(find(model.c));
    mue = findRxnIDs(model, biomassRxnName);

    single_valve_id = findRxnIDs(model, single_valve_reaction)
    glucose = findRxnIDs(model,{glucose});
    o2in = findRxnIDs(model,{'EX_o2(e)'});
    atpm = findRxnIDs(model,{'ATPM'});
    prodind = findRxnIDs(model, target);
    
    if ~isempty(no_KO)
        no_KO = findRxnIDs(model,no_KO);
    else
        no_KO = []
    end
    
    exchange_ind = strmatch('EX',model.rxns);
    DM_rxns = strmatch('DM_',model.rxns);
    exchange_ind = union(exchange_ind,DM_rxns);
    exchange_ind = setdiff(exchange_ind,single_valve_id);

    % Will only find reactions starting with transport ~648 in iJO
    % transpRxns = strmatch('Transport',model.subSystems);
    % Will find all reactions with transport in name ~760 in iJO
    A = regexp(model.subSystems,regexptranslate('wildcard','*Transport*'),'once');
    matched = [];
    for i = 1:length(A)
        if A{i} == 1
            matched(i) = 1;
        else
            matched(i) = 0;
        end
    end
    transpRxns = find(matched)';

    cna_model = CNAcobra2cna(model);
    cna_model_rxn_names = model.rxns;
               
    noKO = [transpRxns; mue; prodind; exchange_ind; no_KO]';

    % remove o2
    noKO = setdiff(noKO,o2in);
    
    if ~isempty(find(noKO==0))
        total_solution.error = 0;
        total_solution.error_message = 'One or more of the reactions was not found';
        disp('One or more of the reactions was not found');
        return;
    end

    st= cna_model.stoichMat; % no external metabolites in cna_model
    irr= cna_model.reacMin >= 0;
    reac_off= union(reac_off, del_exchanges);
    exchange_ind = setdiff(exchange_ind, del_exchanges);
    reac_off = setdiff(reac_off,prodind);   % ensure that the product is not blocked

    product_index = find(model.S(:,findRxnIDs(model,target)));
    if length(product_index) > 1
        total_solution.error = 0;
        total_solution.error_message = 'More than one product index found for the target exchange';
        disp('More than one product index found for the target exchange');
        return;
    end

    glucose_uptake_limit = -model.lb(glucose);
    min_atpm= model.lb(atpm);
    num_seeds= 1;
    separate_reac= unique([noKO, atpm, glucose, mue, o2in]);
    separate_reac(end+1)= 0; % placeholder for index of the export reaction


    %% Begin solution
    %% either initialize or load variables
    q= size(st, 1);
    max_prod= NaN(1, q);
    max_mue= NaN(1, q);
    max_prod_rd= NaN(1, q);
    max_mue_rd= NaN(1, q);
    comp_time= zeros(1, q);
    infeasible= false(1, q);
    is_cs= false(1, q);
    is_des= false(1, q);
    proper_mcs= cell(1, q);
    mcs= cell(1, q);
    is_cs2= false(1, q);
    is_min2= false(1, q);
    is_des2= false(1, q);
    full_sol= cell(1, q);
    gluc_up= NaN(1, q);
    rd_rat= cell(1, q);
    irrev_rd_rat= cell(1, q);
    sub_rat= cell(1, q);
    val_red= NaN(1, q);
    fvalb= cell(1, q);
    fvaub= cell(1, q);
    flux_lb= cell(1, q);
    flux_ub= cell(1, q);

    production= find(st(product_index, exchange_ind));

    if isempty(production)
        stex= [st, zeros(size(st, 1), 1)]; % export reaction is the last reaction
        irrex= [irr; true];
        stex(product_index, end)= -1;
        production= size(st, 2) + 1;
    else
        stex= st;
        irrex= irr;
        production= exchange_ind(production);
    end

    separate_reac(end)= production;

    q= size(stex, 2);
    cgp= Cplex();
    cgp.Param.emphasis.numerical.Cur= 1;
    cgp.Param.preprocessing.presolve.Cur= 0;
    cgp.Param.simplex.tolerances.optimality.Cur= 1e-9;
    cgp.Param.simplex.tolerances.feasibility.Cur= 1e-9;
    cgp.Model.A= stex;
    cgp.Model.ub= Inf(q, 1);
    cgp.Model.lb= -Inf(q, 1);
    cgp.Model.lb(irrex)= 0;
    cgp.Model.lhs= zeros(size(st, 1), 1);
    cgp.Model.rhs= zeros(size(st, 1), 1);
    cgp.Model.lb(glucose)= -glucose_uptake_limit;
    cgp.Model.lb(atpm)= min_atpm;
    cgp.Model.lb(reac_off)= 0;
    cgp.Model.ub(reac_off)= 0;
    cgp.DisplayFunc= [];
    cgp.Model.sense= 'maximize';
    cgp.Model.obj= zeros(q, 1);
    cgp.Model.obj(production)= 1;
    res= cgp.solve();
    if res.status == 1
        max_prod(product_index)= res.objval;
        gluc_up(product_index)= -res.x(glucose);
    elseif res.status == 2
        max_prod(product_index)= Inf;
        gluc_up(product_index)= -res.x(glucose);
    else
        total_solution.error = 1;
        total_solution.error_message = 'Unexpected solution status';
        disp('Unexpected solution status');
        return;
    end
    cgp.Model.obj= zeros(q, 1);
    cgp.Model.obj(mue)= 1;
    res= cgp.solve();
    if res.status == 1
        max_mue(product_index)= res.objval;
    else
        total_solution.error = 1;
        total_solution.error_message = 'Unexpected solution status';
        disp('Unexpected solution status');
        return;
    end

    if isinf(max_prod(product_index)) || (max_prod(product_index) < 1e-9)
        total_solution.error = 1;
        total_solution.error_message = 'metabolite product_index is not limited by glucose or not producible from glucose';
        disp('metabolite is not limited by glucose or not producible from glucose');
        return;
    end

    prod_min_yield= max_prod(product_index)/gluc_up(product_index)*production_stage_min_yield;
    if isempty(sub_rat{product_index})
        des= zeros(0, size(stex, 2));
        des(1, [mue, glucose])= [-1, -production_stage_min_mue]; % growth yield
        des(2, glucose)= -1;
        des(3, [production, glucose])= [-1, -prod_min_yield];
        des(4, atpm)= -1;
        db= [0; glucose_uptake_limit; 0; -min_atpm]; % growth yield

        coli_full= struct();
        coli_full.stoichMat= stex;
        coli_full= CNAgenerateMFNetwork(coli_full);
        coli_full.reacMin(~irrex)= -Inf;
        coli_full.reacMin(reac_off)= 0;
        coli_full.reacMax(reac_off)= 0;
        lpub= coli_full.reacMax;
        lplb= coli_full.reacMin;
        %     [ind_i, ind_j, val]= find(stex);
        %     % using emphasizeFeasibility here yields inaccurate results so
        %     % emphasizeAccuracy is used instead
        %     fvaj= CplexFVA.fva_mt(size(stex, 1), ind_i-1, ind_j-1, val, lplb, lpub, emphasizeAccuracy);
        lhs= [zeros(coli_full.nums, 1); -Inf(length(db), 1)];
        rhs= [zeros(coli_full.nums, 1); db];
        [ind_i, ind_j, val]= find([stex; des]);
        % do FVA with des so that the results can be used for validation later
        fvaj= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);

        if fvaj(3)
            total_solution.error = 1;
            total_solution.error_message = 'FVA infeasible.';
            disp('FVA infeasible.');
            infeasible(product_index)= true;
            return;
        end

        clear ind_i ind_j val;
        fvalb{product_index}= fvaj(1);
        fvaub{product_index}= fvaj(2);
        fva_tol= 1e-9;
        fvalb{product_index}(abs(fvalb{product_index}) < fva_tol)= 0;
        fvaub{product_index}(abs(fvaub{product_index}) < fva_tol)= 0;
        same= fvalb{product_index} == fvaub{product_index};
        blocked_fva= same & (fvalb{product_index} == 0);
        separate_reac(blocked_fva(separate_reac))= []; % do not keep blocked reactions separately
        [rd, irrev_rd, sub]= CNAcompressMFNetwork(coli_full, separate_reac, [], 1, 0, 1, find(blocked_fva), 1);
        sub= sparse(sub');
        %     [rd, sub, irrev_rd]= subsets_reduction(stex, irrex, blocked_fva, separate_reac);
        rd= sparse(rd);
        rd_rat{product_index}= rd;
        irrev_rd= logical(irrev_rd);
        irrev_rd_rat{product_index}= irrev_rd;
        sub_rat{product_index}= sub;
    else
        rd= rd_rat{product_index};
        irrev_rd= irrev_rd_rat{product_index};
        sub= sub_rat{product_index};
    end
    muer= find(sub(:, mue));
    glucr= find(sub(:, glucose));
    prodr= find(sub(:, production));
    atpmr= find(sub(:, atpm));

    cgp= Cplex();
    q= size(rd, 2);
    cgp.Param.emphasis.numerical.Cur= 1;
    cgp.Param.advance.Cur= 0;
    cgp.Param.preprocessing.presolve.Cur= 0;
    cgp.Param.simplex.tolerances.optimality.Cur= 1e-9;
    cgp.Param.simplex.tolerances.feasibility.Cur= 1e-9;
    cgp.Model.A= rd;
    cgp.Model.ub= Inf(q, 1);
    cgp.Model.lb= -Inf(q, 1);
    cgp.Model.lb(irrev_rd ~= 0)= 0;
    % here glucose consumption has a negative sign in the reduced system too
    cgp.Model.lb(glucr)= -glucose_uptake_limit;
    cgp.Model.lb(atpmr)= min_atpm;
    cgp.Model.lhs= zeros(size(rd, 1), 1);
    cgp.Model.rhs= zeros(size(rd, 1), 1);
    cgp.Model.sense= 'maximize';
    cgp.DisplayFunc= [];
    cgp.Model.obj= zeros(q, 1);
    cgp.Model.obj(muer)= 1;
    res= cgp.solve();

    if res.status == 1
        max_mue_rd(product_index)= res.objval;
    else
        total_solution.error = 2;
        total_solution.error_message = 'LP_problem';
        disp('LP problem');
        return;
    end
    cgp.Model.obj= zeros(q, 1);
    cgp.Model.obj(prodr)= 1;
    res= cgp.solve();
    if res.status == 1
        max_prod_rd(product_index)= res.objval;
        %     gluc_up(product_index)= -res.x(glucr);
    else
        total_solution.error = 2;
        total_solution.error_message = 'LP_problem';
        disp('LP problem');
        return;
    end

    inh= zeros(0, size(rd, 2));
    inh(1, glucr)= -1;
    inh(2, [prodr, glucr])= [1, prod_min_yield]; % +prod_min_yield because glucose influx has a negative sign!
    inh(3, atpmr)= -1;
    ub= [glucose_uptake_limit; 0; -min_atpm];
    des= zeros(0, size(rd, 2));
    des(1, [muer, glucr])= [-1, -production_stage_min_mue]; % growth yield
    des(2, glucr)= -1;
    des(3, [prodr, glucr])= [-1, -prod_min_yield];
    des(4, atpmr)= -1;
    db= [0; glucose_uptake_limit; 0; -min_atpm]; % growth yield

%         des2= zeros(0, size(rd, 2));
%         des2(1, muer)= -1;
%         db2= -max_mue_rd(product_index)*growth_stage_min_mue; %valve

    des2= zeros(0, size(rd, 2));
    des2(1, muer)= -1;
    des2(2, glucr)= -1;
    des2(3, atpmr)= -1;
    db2= [-max_mue_rd(product_index)*growth_stage_min_mue; glucose_uptake_limit; -min_atpm]; %valve

    [m, n]= size(rd);
    lplb= -Inf(n, 1);
    lplb(irrev_rd ~= 0)= 0;
    lpub= Inf(n, 1);
    lhs= [zeros(m, 1); -Inf(length(db), 1)];
    rhs= [zeros(m, 1); db];
    [ind_i, ind_j, val]= find([rd; des]);
    valj= CplexValidateMCS.validate(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, 0:n-1, 0:n-1, emphasizeAccuracy);
    fvajrd= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
    clear ind_i ind_j val;
    if fvajrd(3) % FVA was aborted
        total_solution.error = 1;
        total_solution.error_message = 'FVA infeasible.';
        disp('FVA infeasible.');
        infeasible(product_index)= true;
        mcs{product_index}= [];
        return;
    end
    if isnan(val_red(product_index))
        val_red(product_index)= max(max(abs(validate_reduction_fva_bounds(sub, fvaj(1), fvaj(2), fvajrd(1),fvajrd(2)))));
    end
    lhs= [zeros(m, 1); -Inf(length(db2), 1)];
    rhs= [zeros(m, 1); db2];
    [ind_i, ind_j, val]= find([rd; des2]);
    fvajrd2= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
    clear ind_i ind_j val;
    if fvajrd2(3) % FVA was aborted
        total_solution.error = 1;
        total_solution.error_message = 'FVA infeasible.';
        disp('FVA infeasible.');
        infeasible(product_index)= true;
        mcs{product_index}= [];
        return;
    end
    flux_lb2= fvajrd2(1);
    flux_ub2= fvajrd2(2);
    flux_lb2(isinf(flux_lb2))= -2000;
    flux_ub2(isinf(flux_ub2))= 2000;

    flux_lb{product_index}= fvajrd(1);
    flux_ub{product_index}= fvajrd(2);
    flux_lb{product_index}(isinf(flux_lb{product_index}))= -2000;
    flux_ub{product_index}(isinf(flux_ub{product_index}))= 2000;
    cuts= valj(1)' == 0;
    ind= find(cuts == true);
    ind= ind(flux_lb{product_index}(ind) > 0);
    if ~isempty(ind)
        if any(flux_lb{product_index}(ind) > 1e-9)
            total_solution.error = 1;
            total_solution.error_message = 'ERROR: Encountered non-essential reaction with lower bound > 1e-9; skipping.';
            disp('ERROR: Encountered non-essential reaction with lower bound > 1e-9; skipping.');
            return;
        end
        total_solution.warning = 1;
        total_solution.error_message = 'Warning: Encountered non-essential reaction with positive lower bound <= 1e-9; setting bound to 0.';
        disp('Warning: Encountered non-essential reaction with positive lower bound <= 1e-9; setting bound to 0.');
        flux_lb{product_index}(ind)= 0;
    end
    cuts(any(sub(:, noKO), 2))= false;
    
    % disp(single_valve_reaction)
    % if ~strcmp(single_valve_reaction,'')
    %     no_valve = setdiff(1:n, find(sub(:, single_valve_id)));
    % else
    %     no_valve = [];
    % end
    if ~strcmp(single_valve_reaction,'')
        no_valve = setdiff(1:n, find(sub(:, single_valve_id)));
    elseif ~isempty(no_valve2)
        [no_valve,cols,vals] = find(sub(:, no_valve2));
    else
        no_valve = [];
    end

    no_valve = unique(no_valve);

    s= 1;
    while s <= num_seeds %&& (~is_cs(product_index) || ~is_des(product_index)) && ~infeasible(product_index)
        %     obj= ConstrainedMinimalCutSetsEnumerator(rd, irrev_rd, [], inh, ub, cuts, flux_lb{product_index}, flux_ub{product_index}, des, db);
        %     obj.cpx.setParam(cplex_inner.DoubleParam.EpGap, 0.98);

        obj= Orthogonal_cMCSEnumerator(rd, irrev_rd, [], inh, ub, cuts,...
            flux_lb{product_index}, flux_ub{product_index}, des, db, flux_lb2, flux_ub2, des2, db2, no_valve);

        if distributed
            % obj.cpx.readCopyVMConfig(vmconfig); % set virtual machine configuration (or read from a file with obj.cpx.readVMConfig(vmcfile);)
            obj.cpx.copyVMConfig(vmconfig); % set virtual machine configuration (or read from a file with obj.cpx.readVMConfig(vmcfile);)
        end
        
        % Set-up the problem specifics
        if strcmp(objective_type,'unweighted')
            %     This is the case 1 where we want only one valve
            if strcmp(valve_constraint,'equal')
                % Equality constraint
                obj.cpx.addEq(obj.diff_expr, number_of_valves); %implements Sum(yold)-sum(ynew)=1
            elseif strcmp(valve_constraint,'less_than_or_equal')
                %Inequality constraint
                obj.cpx.addLe(obj.diff_expr, number_of_valves); %Implementing const. Sum(yold)-sum(ynew)<=1
            else
                total_solution.error = 3;
                total_solution.error_message = 'Unknown valve_constraint';
                error('Unknown valve_constraint');
            end

            % unweighted objective
            obj= set_objective_function(obj); %without weight
        elseif strcmp(objective_type,'weighted')
            %This is the case 2 where we can have more valves by adjusting weight
            obj= set_objective_function(obj,objective_weight);
        else
            total_solution.error = 3;
            total_solution.error_message = 'Unrecognized objective type';
            error('Unrecognized objective type');
        end

        % calculate cutsets
        % obj.cpx.exportModel('xyz.sav')
        [obj, mcs{product_index}, status, ct]= shortest_minimal_cut_sets(obj, n, num_evs_lim, timeout, solution_type, 0, MILPparameters);
        total_solution.status = status;
        % %Following        
        comp_time(product_index)= comp_time(product_index) + ct;
        infeasible(product_index)= status.equals(cplex_inner.Status.Infeasible);

        if ~isempty(mcs{product_index}) && ~infeasible(product_index)
            [is_cs(product_index), ~, is_des(product_index)]= validate_mcs2(rd, irrev_rd, inh, ub, mcs{product_index},...
                des, db, lplb, lpub);
            if is_cs(product_index) && is_des(product_index)
                full_sol{product_index}= obj.cpx.getValues(obj.mipmat);
                proper_mcs{product_index}= make_minimal_cs(rd, irrev_rd_rat{product_index}, inh, ub, mcs{product_index}, lplb, lpub);
                [is_cs2(product_index), is_min2(product_index), is_des2(product_index)]= validate_mcs2(rd, irrev_rd_rat{product_index},...
                    inh, ub, proper_mcs{product_index}, des, db, lplb, lpub);
            end
        end


        s= s + 1;
    end


    pause(1) % Ensure that this output comes after the Java stdout

    %% check if cut set is minimal (fulfilled when sum(mcs{product_index}) == sum(proper_mcs{product_index}))
    a = sum(mcs{product_index});
    b = sum(proper_mcs{product_index});
    if a > b
        disp('MCS not equal to proper MCS')
    end

    try
        fv2on= round(obj.cpx.getValues(obj.fv2_on));
    catch exception
        total_solution.error = 4;
        total_solution.error_message = 'Failed to calculate solution';
        disp('Failed to calculate solution')
        disp(target)
        disp(exception)
        disp(target);
        disp(exception.message);
        total_solution.name = target;
        total_solution.solution = 0;
        return
    end
    
    % Print the generic solution (reduced model)
    sys.sub= sub_rat{product_index};
    sys= add_subset_names(sys, cna_model_rxn_names);
    char(strcat(num2str(fv2on(mcs{product_index})), ':', sys.sub_name(mcs{product_index})))

    % Loop through the alternative solutions (expanded reduced model)
    cmcsex= expand_mcs(mcs{product_index}, sub_rat{product_index});
    for j= 1:size(cmcsex, 2)
        ind= find(cmcsex(:, j));
        ko_indices{j} = ind;
        strcat(num2str(ind), '|', cna_model_rxn_names(ind))   % This ends up spamming for many results. Supressing output
    end

    comp_growth_stage_cutset = [];
    comp_production_stage_cutset = [];
    valves = [];
    cutset_indices = find(mcs{product_index});
    mcs_size = size(cutset_indices,1);
    for mcs_index = 1:mcs_size
        cut = cutset_indices(mcs_index);
        if fv2on(cut) == 1
            comp_production_stage_cutset = [comp_production_stage_cutset cut];
            valves = [valves cut];
        elseif fv2on(cut) == 0
            comp_production_stage_cutset = [comp_production_stage_cutset cut];
            comp_growth_stage_cutset = [comp_growth_stage_cutset cut];
        else
            total_solution.error = 9;
            total_solution.error_message = 'Unknown fv2on value';
            disp(fv2on(cut))
            disp('Unknown fv2on value');
            return;
        end
    end

    % Expand the growth mcs
    count = 1;
    for k = 1:length(comp_growth_stage_cutset)
        temp_indices = find(sys.sub(comp_growth_stage_cutset(k),:));
        growth_mcs_indices(count) = temp_indices(1);
        count = count + 1;
    end

    % In the case where all the knockouts are valves, set an empty matrix since it is not set above
    if count == 1
        growth_mcs_indices = []
    end

    % Expand the production mcs
    count = 1;
    for k = 1:length(comp_production_stage_cutset)
        temp_indices = find(sys.sub(comp_production_stage_cutset(k),:));
        production_mcs_indices(count) = temp_indices(1);
        count = count + 1;
    end

    % Expand the valve mcs
    count = 1;
    for k = 1:length(valves)
        temp_indices = find(sys.sub(valves(k),:));
        valve_indices(count) = temp_indices(1);
        count = count + 1;
    end

    % Print the result
    disp('Total knockout set');
    printRxnFormula(model,model.rxns(production_mcs_indices));
    disp('')
    disp('Valves');
    printRxnFormula(model,model.rxns(valve_indices));

    total_solution.name = target;
    total_solution.growth_mcs_indices = growth_mcs_indices;
    total_solution.production_mcs_indices = production_mcs_indices;
    total_solution.valve_indices = valve_indices;
    total_solution.solution = 1;
    
    if plot_envelope
        % Plot the three envelopes: wt, growth, production
        f = figure;
        hold on;

        % The substrate uptake must be fixed to ensure correct envelopes
        % Alternatively, we can calculate the yield explicitely
        model = changeRxnBounds(model,model.rxns(glucose),model.lb(glucose),'u');
        if isempty(del_exchanges)
            model_temp = model;
        else            
            model_temp = changeRxnBounds(model,model.rxns(del_exchanges),0,'b');
        end
        
        productionEnvelope(model_temp,model.rxns(growth_mcs_indices),'r',target,biomassRxnName);
        productionEnvelope(model_temp,model.rxns(production_mcs_indices),'blue',target,biomassRxnName);
        productionEnvelope(model_temp,{},'k',target,biomassRxnName);
        legend('growth stage','production stage','wild-type');
        title(strrep(target,'EX_',''));
    end
    pause(1);
end

