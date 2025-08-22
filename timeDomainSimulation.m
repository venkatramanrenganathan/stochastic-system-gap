%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates to obtain distance between two linear stochastic
% systems in the time domain setting
%
% Copyrights @ 2025
% 
% Authors: Venkatraman Renganathan 
%          Cranfield University, United Kingdom.
%
% Email: v.renganathan@cranfield.ac.uk
%
% Date last updated: 22 August, 2025.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a fresh start
clear; close all; clc;

% set properties for plotting
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
addpath(genpath('src'));

%% Set simulation parameters
n = 2;       % state dimension 
m = 1;       % control dimension
l = 1;       % output dimension
d = 4;       % number of uncertain parameters
N = 50;      % number of samples 

% Define matrices for nominal system models
% System 1
A1_0 = [0 1; -2 -0.5];
B1_0 = [0; 1];
C1_0 = [1 0];
D1_0 = 0;
% System 2
A2_0 = [-3.2178  1.2354; -1.7812   -2.6507];
B2_0 = [0; 1];
C2_0 = [1 0];
D2_0 = 0;

% Form the nominal state space models
system1NominalModel = ss(A1_0, B1_0, C1_0, D1_0);
system2NominalModel = ss(A2_0, B2_0, C2_0, D2_0);

% Nominal gap metric
gapBetweenNominalModels = gapmetric(system1NominalModel, system2NominalModel);

% Initialize the random number generator to make the results repeatable
rng(0,'twister');

% Set affine perturbation directions
for k = 1:d
    A1_k(:,:,k) = randn(n);
    B1_k(:,:,k) = randn(n,m);
    C1_k(:,:,k) = randn(l,n);

    A2_k(:,:,k) = randn(n);
    B2_k(:,:,k) = randn(n,m);
    C2_k(:,:,k) = randn(l,n);
end

%% Random Gaussian parameters
meanThetaSystem1 = 0.01;
meanThetaSystem2 = 0.05;
sigmaThetaSystem1 = 0.01;
sigmaThetaSystem2 = 0.05;
thetaSystem1 = sigmaThetaSystem1.*randn(d, N) + meanThetaSystem1;
thetaSystem2 = sigmaThetaSystem2.*randn(d, N) + meanThetaSystem2;

% Storage for N samples of perturbed state-space models for each system
system1PerturbedModels = cell(1, N);
system2PerturbedModels = cell(1, N);

% Generate N samples of perturbed state-space models for each system
for i = 1:N

    % Nominal Matrices
    A1 = A1_0; B1 = B1_0; C1 = C1_0;
    A2 = A2_0; B2 = B2_0; C2 = C2_0;

    % Add affine perturbations to nominal system using randomized theta
    for k = 1:d
        A1 = A1 + thetaSystem1(k,i)*A1_k(:,:,k);
        B1 = B1 + thetaSystem1(k,i)*B1_k(:,:,k);
        C1 = C1 + thetaSystem1(k,i)*C1_k(:,:,k);
        A2 = A2 + thetaSystem2(k,i)*A2_k(:,:,k);
        B2 = B2 + thetaSystem2(k,i)*B2_k(:,:,k);
        C2 = C2 + thetaSystem2(k,i)*C2_k(:,:,k);
    end

    % Store the random perturbed model
    system1PerturbedModels{i} = ss(A1, B1, C1, D1_0);
    system2PerturbedModels{i} = ss(A2, B2, C2, D2_0);
end

%% Compute cost matrix using gap metric
fprintf('Computing %d x %d gap metric based cost matrix...\n', N, N);
gapCost = zeros(N, N);
for i = 1:N
    for k = 1:N
        try
            gapCost(i,k) = gapmetric(system1PerturbedModels{i}, system2PerturbedModels{k});
        catch
            % Set max gap = 1 if system is unstable or error
            gapCost(i,k) = 1; 
        end
    end
end

% Setup the optimal transport linear program
% Placeholder for constraints
Aeq = zeros(2*N, N^2);          
beq = ones(2*N,1)/N;

% Constraints 
for i = 1:N
    Aeq(i, (i-1)*N + (1:N)) = 1;
end

for k = 1:N
    Aeq(N+k, k:N:end) = 1;
end

% Set the objective function
objectiveFunction = gapCost(:);

% Solve the optimal transport linear program
options = optimoptions('linprog','Display','none');
[xopt, fval] = linprog(objectiveFunction, [], [], Aeq, beq, zeros(N^2,1), [], options);

% Get the evaluated gap between systems
type1WassersteinSystemGap = fval;

% Compute upper bound using support distance
systemsGapUpperBound = max(gapCost(:));

% Compute the mean deviation from nominal model for each system
system1DeviationFromNominal = mean(arrayfun(@(i) gapmetric(system1PerturbedModels{i}, system1NominalModel), 1:N));
system2DeviationFromNominal = mean(arrayfun(@(j) gapmetric(system2PerturbedModels{j}, system2NominalModel), 1:N));
systemsGapLowerBound = max(gapBetweenNominalModels - system1DeviationFromNominal - system2DeviationFromNominal, 0);

%% Display Results
fprintf('Distance Lower Bound = %.4f\n', systemsGapLowerBound);
fprintf('Type-1 Wasserstein Distance = %.4f\n', type1WassersteinSystemGap);
fprintf('Distance Upper Bound = %.4f\n', systemsGapUpperBound);

% %% Plot
% figure;
% bar([systemsGapLowerBound, type1WassersteinSystemGap, systemsGapUpperBound], 'FaceColor', 'flat');
% xticklabels({'Lower Bound', '$\mathrm{dist}_{\Sigma_1, \Sigma_2, \delta_{g}}$', 'Upper Bound'});
% ylabel('Distance Between $\Sigma_1$ and $\Sigma_2$'); 
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 8);
% set(a, 'linewidth', 8);
% set(a, 'FontSize', 45);
% set(gca,'fontweight','bold');