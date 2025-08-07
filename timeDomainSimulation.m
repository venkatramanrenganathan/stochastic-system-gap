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
% Date last updated: 06 August, 2025.
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

%% Parameters
n = 2; m = 1; l = 1;
p = 4;        % number of uncertain parameters
N = 50;      % number of samples 

%% Nominal systems (distinct)
A1_0 = [0 1; -2 -0.5];
B1_0 = [0; 1];
C1_0 = [1 0];
D1_0 = 0;

A2_0 = [0 1; -1.5 -0.2];
B2_0 = [0; 1];
C2_0 = [1 0];
D2_0 = 0;

%% Affine perturbation directions
rng(42);
for j = 1:p
    A1_j(:,:,j) = randn(n)*0.01;
    B1_j(:,:,j) = randn(n,m)*0.01;
    C1_j(:,:,j) = randn(l,n)*0.01;

    A2_j(:,:,j) = randn(n)*0.01;
    B2_j(:,:,j) = randn(n,m)*0.01;
    C2_j(:,:,j) = randn(l,n)*0.01;
end

%% Random Gaussian parameters
theta1 = randn(p, N);
theta2 = randn(p, N);

%% Generate perturbed state-space models
sys1_set = cell(1, N);
sys2_set = cell(1, N);

for i = 1:N
    A1 = A1_0; B1 = B1_0; C1 = C1_0;
    A2 = A2_0; B2 = B2_0; C2 = C2_0;

    for j = 1:p
        A1 = A1 + theta1(j,i)*A1_j(:,:,j);
        B1 = B1 + theta1(j,i)*B1_j(:,:,j);
        C1 = C1 + theta1(j,i)*C1_j(:,:,j);

        A2 = A2 + theta2(j,i)*A2_j(:,:,j);
        B2 = B2 + theta2(j,i)*B2_j(:,:,j);
        C2 = C2 + theta2(j,i)*C2_j(:,:,j);
    end

    sys1_set{i} = ss(A1, B1, C1, D1_0);
    sys2_set{i} = ss(A2, B2, C2, D2_0);
end

%% Compute cost matrix using gap metric
fprintf('Computing %d x %d gap metric cost matrix...\n', N, N);
C = zeros(N, N);
for i = 1:N
    for j = 1:N
        try
            C(i,j) = gapmetric(sys1_set{i}, sys2_set{j});
        catch
            C(i,j) = 1; % Max gap if unstable or error
        end
    end
end

%% Setup and solve optimal transport LP
f = C(:);
Aeq = zeros(2*N, N^2);
beq = ones(2*N,1)/N;

for i = 1:N
    Aeq(i, (i-1)*N + (1:N)) = 1;
end

for j = 1:N
    Aeq(N+j, j:N:end) = 1;
end

options = optimoptions('linprog','Display','none');
[xopt, fval] = linprog(f, [], [], Aeq, beq, zeros(N^2,1), [], options);
wass_gap = fval;

%% Support distance
supp_gap = max(C(:));

%% Nominal gap metric
sys_nom1 = ss(A1_0, B1_0, C1_0, D1_0);
sys_nom2 = ss(A2_0, B2_0, C2_0, D2_0);
gap_nom = gapmetric(sys_nom1, sys_nom2);

% Mean deviation from nominal
dev1 = mean(arrayfun(@(i) gapmetric(sys1_set{i}, sys_nom1), 1:N));
dev2 = mean(arrayfun(@(j) gapmetric(sys2_set{j}, sys_nom2), 1:N));
gap_lower = abs(gap_nom - dev1 - dev2);

%% Plot
figure;
bar([gap_lower, wass_gap, supp_gap], 'FaceColor', 'flat');
xticklabels({'Lower Bound', '$\mathrm{dist}_{\Sigma_1, \Sigma_2, \delta_{g}}$', 'Upper Bound'});
ylabel('Distance Between $\Sigma_1$ and $\Sigma_2$'); 
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 8);
set(a, 'linewidth', 8);
set(a, 'FontSize', 45);
set(gca,'fontweight','bold');