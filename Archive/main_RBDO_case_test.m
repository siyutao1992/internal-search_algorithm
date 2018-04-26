%% ME 495: HW2: Siyu TAO

function main_RBDO_case_test

%% Workspace Initialization
clearvars; close all;

SEED = 100;
rng(SEED);

% constants def
P = 150*1e3;    % unit: N
T = 2.5;        % unit: mm
B_mu = 750;        % unit: mm
B_sig = 50;        % unit: mm

E = 210000;     % unit: N/mm^2
smax = 400;     % unit: N/mm^2

% design var def
% x1 -- unit: mm
% x2 -- unit: mm

% design range def (for [x1, x2])
lb = [20, 200];
ub = [80, 1000];

% normal sampling setting
N0 = 1e5;   % initial MCS sampling number
N_rho = 1.2;    % incremental ratio
N_sup = 1e5;    % upper limit of sampling number
MC_conv_rel_tol = 1e-3;   % convergence error fraction


% plot settings
res_n = 51;     % number of pts per axis
linewidth = 1.5;    % line width
markersize = 15;    % marker size

% optimization settings
n = 20;      % number of optimization starting points

% item below for ISO ------------
optsettings = struct(...
    'max_front_step', 0.3,...
    'rho_front_step', 0.5,...
    'min_front_step', 1e-4,...
    'max_later_step', 0.3,...
    'rho_later_step', 0.5,...
    'min_later_step', 1e-4,...
    'max_iter', 100,...
    'max_func_eval',1000,...
    'obj_0th_tol', 1e-6,...
    'obj_1st_tol', 1e-6,...
    'inp_tol', 1e-6);

glob_rel_tol = 5e-2;   % global optimality relative tolerance



%% RBDO w/ internal-search algorithm

% functions def
V_N_func = @(x,B_N) 2*pi*T*x(1)*sqrt(B_N.^2+x(2)^2);
s_smax_N_func = @(x,B_N) P*sqrt(B_N.^2+x(2)^2)/(2*pi*T*x(1)*x(2))-smax; 
s_scrit_N_func = @(x,B_N) P*sqrt(B_N.^2+x(2)^2)/(2*pi*T*x(1)*x(2))-...
    pi^2*E*(T^2+x(1)^2)/8./(B_N.^2+x(2)^2);

prob_less_0 = @(rand_vec) mean((rand_vec<0));

% finding RBDO solution

% starting points (sobol sequence)
soboltemp = sobolset(2,'Skip',1e3);
soboltemp = scramble(soboltemp,'MatousekAffineOwen');
x0_all_norm = net(soboltemp,n);
% start normalization of sobol sequence --------------------------
lb_tmp = min(x0_all_norm); ub_tmp = max(x0_all_norm);
x0_all_norm = (x0_all_norm-repmat(lb_tmp, n, 1))./...
    repmat(ub_tmp-lb_tmp, n, 1);
% end normalization of sobol sequence --------------------------
x0_all = repmat(lb, n, 1)+x0_all_norm.*repmat(ub-lb, n, 1);
x0_all(1,:) = [ub(1)+lb(1), ub(2)+lb(2)]/2;

x_sol_all = zeros(n, 2);
f_sol_all = zeros(n, 1);
ex_fg = zeros(n, 1);
all_out = struct('out_res',cell(n,1),'out_info',cell(n,1));

%
for i = 1:n
    [all_out(i).out_res,all_out(i).out_info] = InterSearchFunc_rect(...
    @(x) V_N_func(x, B_mu), x0_all(i,:), lb, ub, ...
    @(x) nonl_ine_RBDO(x), [], optsettings);
    
    ex_fg(i) = all_out(i).out_info.flag;
    x_sol_all(i,:) = all_out(i).out_res.x_sol;
    f_sol_all(i,:) = all_out(i).out_res.f_sol;
    
    disp([num2str(i),'#',num2str(ex_fg(i,:))]);
end
% filter the converged results
filter = ex_fg>=0;
x0_all = x0_all(filter,:);
x_sol_all = x_sol_all(filter,:); 
f_sol_all = f_sol_all(filter,:);
all_out = all_out(filter);
% sort the results
[f_sol_all, sortI] = sort(f_sol_all);
if abs((f_sol_all(1)-f_sol_all(2))/f_sol_all(1))>glob_rel_tol
    warning('Global optimality checking not passed !');
end
x_sol_all = x_sol_all(sortI,:);
all_out = all_out(sortI);
n_conv = length(f_sol_all);

% calc the n_iter, n_feval, g1 and g2 at all results
n_iter_all = zeros(n_conv, 1);
n_feval_all = zeros(n_conv, 1);
g_opt_all = zeros(n_conv, 2);
for i = 1:n_conv
    [g_opt_all(i,:), ~] = nonl_ine_RBDO(x_sol_all(i,:));
    n_iter_all(i) = all_out(i).out_res.num_iter;
    n_feval_all(i) = all_out(i).out_res.num_func;
end
% assemble all the results to matrix
result_ass_mat = [x0_all, x_sol_all, f_sol_all, g_opt_all, n_iter_all, n_feval_all];

%


% save the results
save(strcat('RBDO_test_case_result',num2str(SEED),'.mat'),...
    'SEED','x0_all', 'x_sol_all', 'g_opt_all','f_sol_all', ...
    'n_iter_all','n_feval_all', 'result_ass_mat');


%%
%% nested functions

    function [c, ceq] = nonl_ine_RBDO(x)
        c = zeros(1,2);
        c(1) = 0.99-QMC_simulator(x,s_smax_N_func, prob_less_0);
        c(2) = 0.99-QMC_simulator(x,s_scrit_N_func, prob_less_0);

        ceq = [];
    end

    function QMC_res = QMC_simulator(x, out_rand_func,formula_func)
        N = N0;
        B_N = qmc(B_mu, B_sig, N);
        if sum(B_N<0) > 0
            warning('Some samples are invalid (B<0) !');
        end
        out_rand = out_rand_func(x,B_N);
        QMC_res = formula_func(out_rand);

        while(true)
            N = round(N*N_rho);
            B_N = qmc(B_mu, B_sig, N);
            MC_res_mark = QMC_res;
            out_rand = out_rand_func(x,B_N);
            QMC_res = formula_func(out_rand);
            if abs((QMC_res-MC_res_mark)/QMC_res) < MC_conv_rel_tol
                break;
            elseif N>N_sup
                %warning('Too many samples required!');
                break;
            end 
        end
    end

end

function B_N = qmc(mu, sig, N) % quasi-MC sampling

qmc=sobolset(1,'Skip',1e3);
qmc = scramble(qmc,'MatousekAffineOwen');
qmc=net(qmc,N);
z_qmc=norminv(qmc,0,1);
B_N = z_qmc*sig+mu;

end



