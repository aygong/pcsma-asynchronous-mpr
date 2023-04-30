function [R, T, parameters, iter] = policy_iteration(mode)
% Apply the policy iteration
% mode: 'upper_bound' or 'heuristic_design'
% 'upper_bound'     : find an upper bound
% 'heuristic_design': find a heuristic design
% Declare global variables
% See main.m
global N gamma c Lambda epsilon max_iter

% Compute the state transition probabilities
% beta : the state transition probabilities
% n1   : the state at slot t
% n2   : the state at slot t + 1
% a    : the action at slot t
% mu   : the probability with which a parameterized policy chooses
%        action a at state n
% p_end: the probability with which every ongoing transmission at slot t
%        is independently completed at the end of slot t
beta = sym(zeros(N+1, N+1));
p_end = 1 / Lambda;
for n1 = 0:N
    for n2 = 0:N
        if n1 < c
            % If a silent user detects n <= c - 1 ongoing transmissions
            % at the beginning of a slot, this user will begin
            % a transmission with probability 0 <= p_n < 1
            for a = max(0, n2-n1):N-n1
                mu = nchoosek(N - n1, a) ...
                    * sym(['p', num2str(n1)])^a * (1 - sym(['p', num2str(n1)]))^(N - n1 - a);
                beta(n1+1, n2+1) = beta(n1+1, n2+1) ...
                    + mu * nchoosek(n1 + a, n1 + a - n2) * p_end^(n1 + a - n2) * (1 - p_end)^n2;
            end
        else
            % Otherwise, this silent user will begin a transmission with probability p_n = 0
            if n1 >= n2
                beta(n1+1, n2+1) = nchoosek(n1, n1 - n2) * p_end^(n1 - n2) * (1 - p_end)^n2;
            end
        end
    end
end

% Compute the modified reward at state n
% rn: the modified reward at state n
% n : the state at slot t
rn = sym(zeros(N+1, 1));
for n = 0:c-1
    for a = 1:gamma-n
        mu = nchoosek(N - n, a) ...
            * sym(['p', num2str(n)])^a * (1 - sym(['p', num2str(n)]))^(N - n - a);
        rn(n+1) = rn(n+1) + Lambda * a * mu;
    end
    if strcmp(mode, 'heuristic_design')
        for a = gamma-n+1:N-n
            mu = nchoosek(N - n, a) ...
                * sym(['p', num2str(n)])^a * (1 - sym(['p', num2str(n)]))^(N - n - a);
            rn(n+1) = rn(n+1) - 2 * Lambda * n * mu;
        end
    end
end

% Initialize two matrices of zeros
% R_heur: the long-term average modified reward
% T_heur: the network throughput
R = zeros(1, max_iter);
T = zeros(1, max_iter);
% Initialize the parameter vector
parameters = zeros(max_iter, c);
parameters(1, 1) = gamma / N;
% Initialize the relative value
v = zeros(N+1, max_iter);

% Run the policy iteration algorithm
for k = 1:max_iter
    % Substitute p^{(k)} into the state transition probabilities
    beta_k = beta;
    for n = 0:c-1
        beta_k = subs(beta_k, sym(['p', num2str(n)]), parameters(k, n+1));
    end
    % Substitute p^{(k)} into the modified reward at state n
    rn_k = rn;
    for n = 0:c-1
        rn_k = subs(rn_k, sym(['p', num2str(n)]), parameters(k, n+1));
    end
    
    % Build the balance equations
    pi_k = sym(zeros(1, N+1));
    for n = 0:N
        pi_k(n+1) = sym(['pi', num2str(n)]);
    end
    one_side = pi_k * beta_k;
    equations = sym(zeros(1, N+1));
    equations(1:N) = one_side(1:N) - pi_k(1:N);
    equations(N+1) = 1 - sum(pi_k);
    % Solve the steady state probabilities
    solutions = solve(equations, pi_k);
    pi_k = double(vpa(struct2cell(solutions)));
    
    % Compute the long-term average modified reward under p^{(k)}
    R(k) = pi_k' * rn_k;
    % Compute the network throughput under p^{(k)}
    T(k) = - throughput_analysis(parameters(k, :));
    
    % Compute the relative value under p^{(k)}
    A = eye(N+1) - beta_k;
    B = rn_k - R(k);
    v(:, k) = pinv(A) * B;
    
    % Update the new parameter vector
    for n = 0:1:c-1
        % Compute the objective function
        objective = rn(n+1) + beta(n+1, :) * v(:, k);
        % Take the derivative of the objective function
        derivative = diff(objective);
        % Find the feasible extremums of the objective function
        solutions = vpasolve(derivative);
        solutions = solutions(imag(solutions) == 0);
        solutions = solutions(and(0 <= solutions, solutions < 1));
        % Add the endpoints of the interval [0, 1]
        solutions = vertcat(solutions, [0; 1]);
        % Update the new parameter
        max_objective = -1e+20;
        for i = 1:length(solutions)
            if subs(objective, solutions(i)) >= max_objective
                max_objective = subs(objective, solutions(i));
                parameters(k+1, n+1) = solutions(i);
            end
        end
    end
    
    % Check the stopping rule
    if sum(abs(parameters(k+1, :) - parameters(k, :))) < epsilon
        fprintf('-> Stop at iteration %d\n', k);
        iter = k;
        break
    end
    
    % Display the number of iterations
    if mod(k, 10) == 0
       fprintf('-> Iteration: %d\n', k);
    end
end