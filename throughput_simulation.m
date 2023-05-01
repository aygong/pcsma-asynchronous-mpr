function R = throughput_simulation(parameters)
% Simulate the network throughput
% Declare global variables
% See main.m
global N gamma c Lambda
global simu_indept simu_slots

% Set the multiple-packet reception channel
channel = zeros(1, N+1);
channel(2:gamma+1) = 1;
% Set the transmission probabilities
p = zeros(1, N+1);
p(1:c) = parameters;

user_index = 1:N;
R = zeros(1, simu_indept);

parfor SI = 1:simu_indept
    % Run independent numerical experiments
    % ongoing  : the status of each user
    %            0 (sensing) or 1 (ongoing)
    % status   : the status of each transmission
    %            0 (unsuccessful) or 1 (successful)
    % length   : the total packet length of each user
    % remaining: the remaining packet length of each user
    % success  : the successful packet length of each user
    ongoing = zeros(1, N);
    status = zeros(1, N);
    lengths = zeros(1, N);
    remaining = zeros(1, N);
    success = zeros(1, N);
    
    for si = 1:simu_slots
        % Simulate consecutive time slots in each experiment
        % Check completed transmissions
        indices = user_index(ongoing == 1);
        remaining(indices) = remaining(indices) - 1;
        indices = intersect(indices, user_index(remaining == 0));
        ongoing(indices) = 0;
        indices = intersect(indices, user_index(status == 1));
        success(indices) = success(indices) + lengths(indices);
        
        % Simulate random access
        if sum(ongoing) < gamma
            access = rand(1, N) < p(sum(ongoing)+1);
            indices = user_index((1 - ongoing) .* access == 1);
            ongoing(indices) = 1;
            status(indices) = 1;
            lengths(indices) = natural_geornd(1 / Lambda, 1, length(indices));
            remaining(indices) = lengths(indices);
        end
        
        % Check unsuccessful transmissions
        if sum(ongoing) > 0 && rand > channel(sum(ongoing)+1)
            status = zeros(1, N);
        end
    end
    
    R(SI) = sum(success) / simu_slots;
end

R = sum(R) / simu_indept;