function [is_cs, is_min, is_des]= validate_mcs2(st, irr, inh, ub, mcs, des, db, lplb, lpub, no_is_min)
% checks with a series of LPs whether the mcs are indeed cut sets, i.e.
% they must make the primal system infeasible, and whether they are
% minimal, i.e. each cut removed from an MCS must result in a feasible
% system
[m, n]= size(st);
num_mcs= size(mcs, 2);
if nargin < 10
  no_is_min= false;
  if nargin < 8 % < 9?
    lplb= -Inf(n, 1);
    lplb(irr ~= 0)= 0;
    lpub= Inf(n, 1);
  end
end
cgp= Cplex();
cgp.Param.emphasis.numerical.Cur= 1;
% use dual simplex and set the feasibility tolerance to its lowest possible
% value to make sure that the constraints which implement the cuts are not violated
% cgp.Param.lpmethod.Cur= 2; % dual simplex
cgp.Param.simplex.tolerances.feasibility.Cur= cgp.Param.simplex.tolerances.feasibility.Min;
cgp.Param.preprocessing.presolve.Cur= 0;
cgp.Model.A= sparse([st; inh]);
cgp.Model.lb= lplb;
cgp.Model.ub= lpub;
cgp.Model.lhs= [zeros(m, 1); -Inf(length(ub), 1)];
cgp.Model.rhs= [zeros(m, 1); ub(:)];
cgp.Model.obj= sparse(1, length(cgp.Model.ub));
cgp.DisplayFunc= [];

if nargin > 6
  cgp2= Cplex();
  cgp2.Param.emphasis.numerical.Cur= 1;
%   cgp2.Param.lpmethod.Cur= 2; % dual simplex
  cgp2.Param.simplex.tolerances.feasibility.Cur= cgp2.Param.simplex.tolerances.feasibility.Min;
  cgp.Param.preprocessing.presolve.Cur= 0;
  cgp2.Model.A= sparse([st; des]);
  cgp2.Model.lb= lplb;
  cgp2.Model.ub= lpub;
  cgp2.Model.lhs= [zeros(m, 1); -Inf(length(db), 1)];
  cgp2.Model.rhs= [zeros(m, 1); db(:)];
  cgp2.Model.obj= sparse(length(cgp2.Model.ub), 1);
  cgp2.DisplayFunc= [];
end

is_cs= false(1, num_mcs);
if no_is_min
  is_min= false(1, num_mcs);
else
  is_min= true(1, num_mcs);
end
is_des= false(1, num_mcs);
for j= 1:num_mcs
  ind= find(mcs(:, j))';
  cgp.Model.lb(ind)= 0;
  cgp.Model.ub(ind)= 0;
  res= cgp.solve();
  is_cs(j)= res.status == 3; % infeasible
  
  if ~no_is_min
    for i= ind % minimality is checked by removing each cut separately
      cgp.Model.lb(i)= lplb(i);
      cgp.Model.ub(i)= lpub(i);
      res= cgp.solve();
      if res.status ~= 1
        is_min(j)= false;
        break;
      end
      cgp.Model.lb(i)= 0;
      cgp.Model.ub(i)= 0;
    end
  end
  
  cgp.Model.lb(ind)= lplb(ind);
  cgp.Model.ub(ind)= lpub(ind);
  
  if nargin > 6
    cgp2.Model.lb(ind)= 0;
    cgp2.Model.ub(ind)= 0;
    res= cgp2.solve();
    is_des(j)= res.status == 1;
    cgp2.Model.lb(ind)= lplb(ind);
    cgp2.Model.ub(ind)= lpub(ind);
  end

  if mod(j, 10) == 0
    fprintf('.');
    if mod(j, 1000) == 0
      fprintf('\n%d', j);
    end
  end
end
if j >= 10
  fprintf('\n');
end
