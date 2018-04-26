function [out_res, out_info] = InterSearchFunc_rect(...
        objfunc, x0, lb, ub, nonl_ine, gradfunc, optsettings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: This function is to normalize (rectify) inputs, but ...
%   the "optsettings" will not be normalized, hence ...
%   please use RELATIVE optimization settings when inputting "optsettings"!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input checking & Initialization
if isempty(gradfunc)
    no_user_grad = true;
else
    no_user_grad = false;
end

if sum((lb-ub)>0)>0
    error('The lb and ub are incompatible!');
end

[be1, d] = size(x0);
if be1 ~= 1
    error('Multiple starting points is not supported yet!');
end

%% normalization
x0_norm = normlz_x(x0);
lb_norm = zeros(1,d);
ub_norm = ones(1,d);

% sample to find the upper and lower limit of obj. It's fine if inaccurate.
if ~isfield(optsettings,'max_func_eval')
    num_func_eval = 1000;
else
    if optsettings.max_func_eval < 500
        error('The allowed maximum number of func eval is too small!');
    end
    num_func_eval = round(0.01*optsettings.max_func_eval);
end
SEED = 123; % This can be a input in the future !
rng(SEED);
soboltemp = sobolset(d,'Skip',SEED*num_func_eval,'Leap',randi(num_func_eval));
soboltemp = scramble(soboltemp,'MatousekAffineOwen');
input_data_norm = net(soboltemp,num_func_eval);
input_data_unnorm = input_data_norm.*repmat(ub-lb,num_func_eval,1)+...
    repmat(lb, num_func_eval,1);
obj_all = zeros(num_func_eval,1);
for i = 1:num_func_eval
    obj_all(i) = objfunc(input_data_unnorm(i,:));
end
obj_lb = min(obj_all);
obj_ub = max(obj_all);

objfunc_norm = @(x_norm) normlz_f(objfunc(unnorm_x(x_norm)));
nonlcon_norm = @(x_norm) nonl_ine(unnorm_x(x_norm));

%% execute the core optimization function
if no_user_grad
    [out_res, out_info] = InterSearchOpt2D(...
        objfunc_norm, x0_norm, lb_norm, ub_norm, ...
        nonlcon_norm, [], optsettings);
else
    gradfunc_norm = @(x_norm) gradfunc(unnorm_x(x_norm))/(obj_ub-obj_lb);

    [out_res, out_info] = InterSearchOpt2D(...
        objfunc_norm, x0_norm, lb_norm, ub_norm, ...
        nonlcon_norm, gradfunc_norm, optsettings);
end

%% unnormalization
out_res.x_sol = unnorm_x(out_res.x_sol);
out_res.f_sol = unnorm_f(out_res.f_sol);
num_iter = length(out_res.f_history);
out_res.x_history = (out_res.x_history).*repmat(ub-lb,num_iter,1)+...
    repmat(lb, num_iter,1);
out_res.f_history = (out_res.f_history).*repmat(obj_ub-obj_lb,num_iter,1)+...
    repmat(obj_lb, num_iter,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% embedded subroutines

    function x_unnorm = unnorm_x(x_norm)
        x_unnorm = x_norm.*(ub-lb)+lb;
    end
    function x_norm = normlz_x(x)
        x_norm = (x-lb)./(ub-lb);
    end
    function f_unnorm = unnorm_f(f_norm)
        f_unnorm = f_norm*(obj_ub-obj_lb)+obj_lb;
    end
    function f_norm = normlz_f(f)
        f_norm = (f-obj_lb)/(obj_ub-obj_lb);
    end

    
end

