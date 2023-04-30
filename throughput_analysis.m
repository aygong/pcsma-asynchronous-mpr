function R = throughput_analysis(parameters)
% Analyze the network throughput
% Declare global variables
% See main.m
global N gamma c Lambda

% Set the transmission probabilities
p = zeros(1, N+1);
p(1:c) = parameters;

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
            for a = max(0, n2 - n1):(N - n1)
                mu = nchoosek(N - n1, a) * p(n1+1)^a * (1 - p(n1+1))^(N - n1 - a);
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

% Build the balance equations
pi = sym(zeros(1, N+1));
for n = 0:N
    pi(n+1) = sym(['pi', num2str(n)]);
end
one_side = pi * beta;
equations = sym(zeros(1, N+1));
equations(1:N) = one_side(1:N) - pi(1:N);
equations(N+1) = 1 - sum(pi);
% Solve the steady state probabilities
solutions = solve(equations, pi);
pi = double(vpa(struct2cell(solutions)));

% Compute the transition probability that there are h' other 
% ongoing transmissions in the next slot when there have been 
% h other ongoing transmissions in the present slot
Xi = zeros(N, N);
for h1 = 0:N-1
    for h2 = 0:N-1
        for j = max(0, h1-h2):h1
            Xi(h1+1, h2+1) = Xi(h1+1, h2+1) ...
                + nchoosek(h1, j) * p_end^j * (1 - p_end)^(h1 - j) ...
                * nchoosek(N - 1 - h1 + j, h2 - h1 + j) ...
                * p(h1 - j + 2)^(h2 - h1 + j) * (1 - p(h1 - j + 2))^(N - 1 - h2);
        end
    end
end

% Initialize the diagonal matrix
mask = zeros(N, N);
mask(1:gamma, 1:gamma) = eye(gamma);

% Compute the reward at state n
% r: the reward at state n
% g: the probability of having h_m other ongoing transmissions in
%    the m-th slot of a given transmission and having fewer than
%    γ other ongoing transmissions in each of the first m − 1 slots
%    of this given transmission, provided that there are h_1 other
%    ongoing transmissions in the first slot
rn = zeros(N+1, 1);
for n = 0:gamma-1
    % Compute the parameterized policy
    mu = zeros(N+1, 1);
    for a = 0:N-n
        mu(a+1) = nchoosek(N - n, a) * p(n+1)^a * (1 - p(n+1))^(N - n - a);
    end
    % Compute the reward that is gained when action a is chosen at state n
    % lambda  : the packet length
    % p_length: the probability that a packet has a length equal to λ slots
    rna = zeros(N+1, 1);
    for a = 1:gamma-n
        % Initialize g
        g = zeros(1, N);
        g(n+a) = 1;
        for lambda = 1:1e+6
            p_length = lambda * p_end * (1 - p_end)^(lambda - 1);
            rna(a+1) = rna(a+1) + p_length * sum(g(1:gamma));
            if p_length < 1e-10
                break
            end
            % Update g
            g = g * mask * Xi;
        end
        rna(a+1) = a * rna(a+1);
    end
    rn(n+1) = dot(mu, rna);
end

% Compute the network throughput
R = - dot(rn, pi);