%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function computeWassersteinLinearProgram will return Wasserstein distance
% given the transort cost using the solution of a linear program
%
% INPUT:
%      transportCost: N x N transport cost matrix
%
% OUTPUT:
%      wassersteinDistance: wasserstein distance in terms of transport cost
%      
%
% Copyrights @ 2025
% 
% Authors: Venkatraman Renganathan 
%          Cranfield University, United Kingdom.
%
% Email: v.renganathan@cranfield.ac.uk
%
% Date last updated: 27 August, 2025.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wassersteinDistance = computeWassersteinLinearProgram(transportCost)
    
    % Solve W1 with uniform marginals via LP: min <D, Gamma>
    % subject to Gamma*1 = 1/N, Gamma^T*1 = 1/N, Gamma >= 0
    [N1,N2] = size(transportCost);
    assert(N1==N2, 'Use equal sample counts on both sides.');
    N = N1;
    objectiveFunction = transportCost(:);

    % Equality constraints
    Aeq = zeros(2*N, N*N);
    beq = ones(2*N,1) / N;
    lb = zeros(N*N,1);

    % Row sums (source)
    for i = 1:N
        idx = (i-1)*N + (1:N);
        Aeq(i, idx) = 1;
    end
    % Column sums (target)
    for j = 1:N
        Aeq(N+j, j:N:end) = 1;
    end

    % Set the options for linprog solver
    opts = optimoptions('linprog','Display','none','Algorithm','dual-simplex');
    
    % Solve the linear program using the linprog solver
    [~, fval, exitflag] = linprog(objectiveFunction, [], [], Aeq, beq, lb, [], opts);
    if exitflag <= 0
        error('linprog failed to solve the linear program (exitflag=%d).', exitflag);
    end

    % Return the Wasserstein-1 cost with given transportCost
    wassersteinDistance = fval; 
end