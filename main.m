clc, clear

%% Declare global variables
% N       : the number of infinitely backlogged users
% gamma   : the multiple-packet reception capability
% c       : the carrier sensing capability
% Lambda  : the average packet length
% epsilon : the threshold for terminating iterations
% max_iter: the maximum number of iterations
global N gamma c Lambda epsilon max_iter
% simu_switch: the simulation switch
% simu_indept: the number of independent numerical experiments
% simu_slots : the number of time slots
global simu_switch simu_indept simu_slots

% Set the network parameters
N = 10;
gamma = 5;
c = 5;
Lambda = 100;

% Set the iteration parameters
epsilon = 1e-10;
max_iter = 100;

% Set the simulation parameters
simu_switch = true;
simu_indept = 8;
simu_slots = 1e+6;

%% Return the network throughput
% Display the network parameters
fprintf('|> N = %d, gamma = %d, c = %d, Lambda = %d\n', N, gamma, c, Lambda);

% Return the upper bound
fprintf('|> Find an upper bound\n');
[R_upp, T_upp, p_upp, iter] = policy_iteration('upper_bound');
if simu_switch
    T_ana = T_upp(iter);
    T_sim = throughput_simulation(p_upp(iter, :));
    results_display(T_ana, T_sim);
end

% Return the network throughput under the heuristic scheme
fprintf('|> Find a heuristic design\n');
[R_heu, T_heu, p_heu, iter] = policy_iteration('heuristic_design');
if simu_switch
    T_ana = T_heu(iter);
    T_sim = throughput_simulation(p_heu(iter, :));
    results_display(T_ana, T_sim);
end

% Return the network throughput under the optimal scheme
fprintf('|> Find an optimal design\n');
ms = MultiStart('FunctionTolerance', 1e-6, 'UseParallel', true, 'Display', 'iter');
opts = optimoptions(@fmincon, 'Display', 'iter-detailed');
gs = GlobalSearch(ms);
problem = createOptimProblem('fmincon', 'x0', p_heu(iter, :), 'objective', ...
    @throughput_analysis, 'lb', zeros(1, c), 'ub', ones(1, c), 'options', opts);
[p_opt, T_opt] = run(gs, problem);
if simu_switch
    T_ana = - T_opt;
    T_sim = throughput_simulation(p_opt);
    results_display(T_ana, T_sim);
end