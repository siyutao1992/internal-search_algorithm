function [out_res, out_info] = InterSearchOpt2D(...
    objfunc, x0, lb, ub, nonl_ine, gradfunc, optsettings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Please normalize ...
%   "objfunc (input & output), x0, lb, ub, nonlcon, gradfunc (if provided)"
%   before using this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=====================================================
%{
    Inputs:
        objfunc:    (fun) objective function
        x0:         (vec) initial starting point
        lb:         (vec) lower search bound
        ub:         (vec) upper search bound
        nonlcon:    (fun) nonlinear constraint function
        gradfunc (optional) :   (fun) gradient function of the objective
                                NOTE: the shape of grad is [1, p] !!
        optsettings: (stt) optimization setting struct
            .max_front_step:    maximum linesearch stepsize forward
            .rho_front_step:    decrease ratio linesearch stepsize forward
            .min_front_step:    minimum linesearch stepsize forward
            .max_later_jump:    maximum lateral jumpsize
            .rho_later_jump:    decrease ratio stepsize lateral
            .min_later_jump:    minimum lateral jumpsize
            .max_iter:              maximum number of iterations
            .max_func_eval:         maximum function evaluations
            .obj_0th_tol:           objective zero-order optimality tolerance
            .obj_1st_tol:           objective 1st-order optimality tolerance
            .inp_tol:               input zero-order optimality tolerance

    Outputs:
        out_res:    (stt) optimization results struct
            num_iter:   (num) number of iterations
            num_func:   (num) number of objective function evaluations
            x_sol:      (vec) optimal input solution
            f_sol:      (num) optimal objective/output
            x_history:  (mat) optimization history for inputs
            f_history:  (vec) optimization history for objective/output
        out_info:   (stt) optimization info summary
            .flag:      convergence reason...
                        "-1" start with infeasible point
                        "0" maximum function evaluation reached
                        "1" convergent off constraint, 0th-order optimality
                        "2" convergent off constraint, 1st-order optimality
                        "3" convergent off constraint, input step optimality
                        "4" convergent at constraint
            .last_front_stepsize:   the last forward stepsize
            .last_later_jumpsize:   the last lateral jumpsize
            .obj_0th_opt:           objective zero-order optimality value
            .obj_1st_opt:           objective 1st-order optimality value
            .inp_opt:               input zero-order optimality value
            
%}
%======================================================
%% Input checking & Initialization
if isempty(gradfunc)
    no_user_grad = true;
else
    no_user_grad = false;
end

if sum((lb-ub)>0)>0
    error('The lb and ub are incompatible!');
end

% dimension extract
[be1, p] = size(x0);
if be1 ~= 1
    error('Multiple starting points is not supported yet!');
end

% optimzation settings regularization
settings = settings_init(optsettings);
max_front_step = settings.max_front_step;
rho_front_step = settings.rho_front_step;
min_front_step = settings.min_front_step;
max_later_step = settings.max_later_jump;
rho_later_step = settings.rho_later_jump;
min_later_step = settings.min_later_jump;
max_iter = settings.max_iter;
max_func_eval = settings.max_func_eval;
obj_0th_tol = settings.obj_0th_tol;
obj_1st_tol = settings.obj_1st_tol;
inp_tol = settings.inp_tol;

% outputs initialization
out_res = struct('num_iter',0,'num_func',0,'x_sol',zeros(1,p),'f_sol',0,...
    'x_history',zeros(max_iter,p), 'f_history',zeros(max_iter,1));
out_info = struct('flag',-1,'last_front_stepsize',max_front_step,...
    'last_later_jumpsize',max_later_step,'obj_0th_opt', 0,...
    'obj_1st_opt', 0, 'inp_opt', 0);

x_history = zeros(max_iter,p);
f_history = zeros(max_iter,1);

% normalization
% THIS SHOULD BE DONE BEFORE USING THIS FUNCTION
% demo
%{
x0_norm = normlz_x(x0);
lb_norm = zeros(1,d);
ub_norm = ones(1,d);
objfunc_norm = @(x_norm) objfunc(unnorm_x(x_norm));
nonlcon_norm = @(x_norm) nonlcon(unnorm_x(x_norm));
%}

%% Iterations
i_iter = 0;
i_func = 0;
x = x0;
f = objfunc(x);
i_func = i_func +1;
while true
    x_prev = x;
    f_prev = f;
    % check starting point feasibility
    if i_iter == 0
        if ~check_feas(x)
            out_info.flag = -1; % This is already set
            break;  % optimization stops due to initial infesaibility
        end
    end
    if i_iter == max_iter || i_func > max_func_eval
        out_info.flag = 0;
        break;  % optimization stops due to max func eval or max iter
    end
        
    
    % forward linesearch
    [front_success, infeas_ahead, grad] = linesearch(x);
    % the "grad" output is useful only if lateral search is executed
    if ~front_success && ~infeas_ahead
        disp('pause');
    end
    later_success = false;
    % lateral jump&linesearch
    if ~front_success && infeas_ahead
        later_dir = [-(null(grad))';(null(grad))'];
        later_step = max_later_step;
        while true
            if later_step < min_later_step
                break;  % lateral search fails, exit lateral search
            end
            for ii = 1:size(later_dir,1)    % search all lateral dirs with this jumpsize
                x_jump = x+later_step*later_dir(ii,:);
                if check_feas(x_jump)
                    i_func = i_func +1;
                    if f-objfunc(x_jump) > 0
                        update_iter(x_jump);
                        record_iter();
                        later_success = true;
                        break;  % lateral search success at jump, exit dirs search
                    end
                    
                    [later_success, ~, ~] = linesearch(x_jump);
                    if later_success
                        break;  % lateral search success at linesearch, exit dirs search
                    end
                end
            end
            if later_success
                break;  % lateral search success, exit lateral search
            else
                later_step = later_step*rho_later_step;
                continue;
            end
        end
    end
    
    % decide if convergent
        
    if front_success && ~infeas_ahead
        if f_prev-f< obj_0th_tol
            out_info.flag = 1;
            break; % convergent off constraint, 0th-order optimality
        elseif norm(grad) < obj_1st_tol
            out_info.flag = 2;
            break; % convergent off constraint, 1st-order optimality
        elseif norm(x-x_prev) < inp_tol
            out_info.flag = 3;
            break; % convergent off constraint, input step optimality
        end
    elseif ~front_success && ~later_success
        out_info.flag = 4;
        break; % convergent at constraint
    end
end    


%% record output
out_res.num_iter = i_iter;
out_res.num_func = i_func;
if i_iter > 0
    out_res.x_sol = x_history(i_iter,:);
    out_res.f_sol = f_history(i_iter,:);
    out_res.x_history = x_history;
    out_res.f_history = f_history;
    %{
    figure;
    hold on;
    for j = 1:i_iter
        plot(x_history(j,1),x_history(j,2),'.','markersize',20);
        text(x_history(j,1),x_history(j,2),num2str(j));
    end
    axis equal
    hold off;
    %}  
end

% "out_info" will be later completed/updated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% embedded subroutines

    function isfeasible = check_feas(x)
        isfeasible = true;
        if sum((x-lb)<0)>0 || sum((ub-x)<0)>0 || sum(nonl_ine(x)>0)>0
            isfeasible = false;
        end
    end
    
    function update_iter(x_in)
        i_iter = i_iter +1;
        x = x_in;
        f = objfunc(x); i_func = i_func +1;
    end

    function record_iter()
        x_history(i_iter,:) = x;
        f_history(i_iter,:) = f;
    end

    function grad = default_grad(x)
        FiniteDiffType = 'centered'; % can be 'forward','backward' or 'center'
        err_tol = 1e-1;
        max_d = 1e-2;
        decre_rho = 0.8;
        min_d = 1e-4;
        
        new_d = max_d;
        grad = one_diff_grad(x, new_d, FiniteDiffType);
        while true
            new_d = new_d * decre_rho;
            if new_d < min_d
                warning('The finite difference grad may be inaccurate!');
                break;
            end
            grad_prev = grad;
            grad = one_diff_grad(x, new_d, FiniteDiffType);
            %disp([grad_prev,grad]);
            %disp(norm(grad-grad_prev));
            if norm(grad-grad_prev) < err_tol
                break;
            end 
        end
        
    end

    function grad = one_diff_grad(x, d, FiniteDiffType)
        grad = zeros(1,p);
        if strcmp(FiniteDiffType, 'forward')
            depart_x = repmat(x,p,1)+d*eye(p);
            obj_cen = objfunc(x);i_func = i_func +1;
            for i = 1:p
                grad(i) = (objfunc(depart_x(i,:))-obj_cen)/d;
                i_func = i_func +1;
            end
        elseif strcmp(FiniteDiffType, 'backward')
            depart_x = repmat(x,p,1)-d*eye(p);
            obj_cen = objfunc(x);i_func = i_func +1;
            for i = 1:p
                grad(i) = (obj_cen-objfunc(depart_x(i,:)))/d;
                i_func = i_func +1;
            end
        elseif strcmp(FiniteDiffType, 'centered')
            forw_depart_x = repmat(x,p,1)+d/2*eye(p);
            back_depart_x = repmat(x,p,1)-d/2*eye(p);
            for i = 1:p
                grad(i) = (objfunc(forw_depart_x(i,:))-...
                    objfunc(back_depart_x(i,:)))/d;
                i_func = i_func +2;
            end
        else
            error('Finite diff type not recognized!');
        end
    end

    function [success, infeas_ahead, grad] = linesearch(x_in)
        
        success = false;
        infeas_ahead = false;
        % determine the gradient direction
        if no_user_grad
            grad = default_grad(x_in);
        else
            grad = gradfunc(x_in);
            if size(grad,1) > 1 
                if size(grad,2) == 1
                    grad = grad';
                else
                    error('The shape of the user''s grad func output is wrong!');
                end
            end
        end
        % 
        front_step = max_front_step;
        while true
            if front_step < min_front_step
                if infeas_ahead == false
                    disp('This is unlikely');
                end
                break;  % line search fails, exit line search
            end
            x_try = x_in-front_step*grad;
            if ~check_feas(x_try)
                infeas_ahead = true;
                front_step = front_step*rho_front_step;
                continue;   % infeasibility ahead, decrease stepsize
            end
            infeas_ahead = false;
            i_func = i_func +1;
            if f-objfunc(x_try) > 0
                update_iter(x_try);
                record_iter();
                success = true;
                break;  % line search success, exit line search
            else
                front_step = front_step*rho_front_step;
                continue;
            end
        end
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% external subroutines
function settings = settings_init(optsettings)

% default values
if ~isfield(optsettings,'max_front_step')
    optsettings.max_front_step = 0.1;
end
if ~isfield(optsettings,'rho_front_step')
    optsettings.rho_front_step = 0.5;
end
if ~isfield(optsettings,'min_front_step')
    optsettings.min_front_step = 1e-4;
end
if ~isfield(optsettings,'max_later_jump')
    optsettings.max_later_jump = 0.1;
end
if ~isfield(optsettings,'rho_later_jump')
    optsettings.rho_later_jump = 0.5;
end
if ~isfield(optsettings,'min_later_jump')
    optsettings.min_later_jump = 1e-4;
end
if ~isfield(optsettings,'max_func_eval')
    optsettings.max_func_eval = 1000;
end
if ~isfield(optsettings,'max_iter')
    optsettings.max_iter = 200;
end
if ~isfield(optsettings,'obj_0th_tol')
    optsettings.obj_0th_tol = 1e-6;
end
if ~isfield(optsettings,'obj_1st_tol')
    optsettings.obj_1st_tol = 1e-6;
end
if ~isfield(optsettings,'inp_tol')
    optsettings.inp_tol = 1e-4;
end

settings = optsettings;
end



