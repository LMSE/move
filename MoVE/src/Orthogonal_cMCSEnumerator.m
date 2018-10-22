classdef Orthogonal_cMCSEnumerator < ConstrainedMinimalCutSetsEnumerator
  properties (SetAccess = 'immutable')
    flux_vars2
    fv2_on
    first_fv2_on;
    diff_expr;
  end

  methods (Access = 'public')
    function obj= Orthogonal_cMCSEnumerator(st, irr, targets, inh, ubi, cuts,...
        flux_lb, flux_ub, des, db, flux_lb2, flux_ub2, des2, db2, no_valve)
      % no_valve is a vector of indices of reactions which are not allowed
      % to be chosen as valves
      if nargin < 15
        no_valve= [];
      end 
      
      obj= obj@ConstrainedMinimalCutSetsEnumerator(st, irr, targets, inh, ubi, cuts,...
        flux_lb, flux_ub, des, db);
      
      [m, n]= size(st);
      numsn= nums(n);
      
      obj.first_fv2_on= obj.mipmat.getNcols(); % Java index
      obj.fv2_on= obj.cpx.boolVarArray(n, strcat('FV2ON', numsn));
      obj.mipmat.addCols(obj.fv2_on);
      
      rind= obj.mipmat.getNrows(); % Java index
      rhs= 2*ones(n, 1);
      rhs(no_valve)= 1;
      obj.mipmat.addRows(ones(n, 1), rhs, [], []);
%       for i= 0:n-1
%         obj.mipmat.setNZs(rind + [i i i], [obj.first_z obj.first_z+n obj.first_fv2_on] + i, [1 1 1]); 
%       end
      obj.mipmat.setNZs(repmat(rind:rind+n-1, 1, 3),...
        [obj.first_z:obj.first_z+2*n-1 obj.first_fv2_on:obj.first_fv2_on+n-1], ones(1, 3*n));
     
      first_fv2= obj.mipmat.getNcols(); % Java index
      obj.flux_vars2= obj.cpx.numVarArray(n, flux_lb2, flux_ub2, strcat('R2_', numsn));
      obj.mipmat.addCols(obj.flux_vars2);
      rind= obj.mipmat.getNrows(); % Java index
      obj.mipmat.addRows([zeros(m, 1); -Inf(length(db2), 1)], [zeros(m, 1); db2], [], []);
      [i, j, val]= find([st; des2]);
      obj.mipmat.setNZs((rind - 1) + i, first_fv2 + int32(j - 1), val);
      
%       for i= 1:n
%         obj.cpx.addLe(obj.flux_vars2(i), obj.cpx.prod(obj.fv2_on(i), flux_ub(i)));
%       end
%       for i= 1:n
%         obj.cpx.addGe(obj.flux_vars2(i), obj.cpx.prod(obj.fv2_on(i), flux_lb(i)));
%       end
      rind= obj.mipmat.getNrows();
      obj.mipmat.addRows([-Inf(n, 1); zeros(n, 1)], [zeros(n, 1); Inf(n, 1)], [], []);
      obj.mipmat.setNZs(rind-1+(1:2*n), [first_fv2:first_fv2+n-1, first_fv2:first_fv2+n-1],...
        ones(2*n, 1));
      obj.mipmat.setNZs(rind-1+(1:2*n), [obj.first_fv2_on:obj.first_fv2_on+n-1, obj.first_fv2_on:obj.first_fv2_on+n-1],...
        -[flux_ub2(:); flux_lb2(:)]);
      
      % this can be used to manually limit the difference between the
      % reaction participations like this:
      % obj.cpx.addLe(obj.diff_expr, 1)
      obj.diff_expr= obj.cpx.diff(obj.cpx.sum(obj.cpx.sum(obj.z_vars),...
        obj.cpx.sum(obj.fv2_on)), n);
        
    end

    function obj= set_objective_function(obj, zv_weight)
      % only sets the objective function if none is present or if the
      % current one is a constant
      if nargin < 2
        obj_expr= obj.cpx.sum(obj.z_vars);
      else
        obj_expr= obj.cpx.sum(obj.cpx.prod(zv_weight, obj.cpx.sum(obj.z_vars)), obj.diff_expr);
      end
      if isempty(obj.cpx.getObjective())
        obj.cpx.addMinimize().setExpr(obj_expr);
      else
        objective= obj.cpx.getObjective();
        if ~objective.getExpr().linearIterator().hasNext()
          objective.setExpr(obj_expr);
        end
      end
    end

  end

end % classdef

function res= nums(n)
res= strtrim(cellstr(num2str((1:n)')));
end
