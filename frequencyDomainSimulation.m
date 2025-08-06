%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates to obtain distance between two linear stochastic
% systems in the frequency domain setting
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


clear; clc;

%% Parameters
N = 500;             % Use N = 1e4 for large scale (memory warning!)
M = 50;              % Frequency points
omega = linspace(0.1, 10, M);

%% Define distinct nominal systems
P1_nom = @(w) 1 ./ (1 + 0.5j * w);
P2_nom = @(w) 1 ./ ((1 + 0.2j * w) .* (1 + 0.7j * w));
P1_vals = P1_nom(omega);
P2_vals = P2_nom(omega);

%% Generate stochastic frequency response samples
perturb_scale = 0.1;
P1_samples = zeros(N, M);
P2_samples = zeros(N, M);

for i = 1:N
    noise1 = perturb_scale * (randn(1, M) + 1j * randn(1, M));
    noise2 = perturb_scale * (randn(1, M) + 1j * randn(1, M));
    P1_samples(i, :) = P1_vals .* (1 + noise1);
    P2_samples(i, :) = P2_vals .* (1 + noise2);
end

%% Inverse stereographic projection
toRiemann = @(z) [real(z) ./ (1 + abs(z).^2); ...
                  imag(z) ./ (1 + abs(z).^2); ...
                  abs(z).^2 ./ (1 + abs(z).^2)];

P1_R = zeros(3, N, M);
P2_R = zeros(3, N, M);
R1_nom = zeros(3, M);
R2_nom = zeros(3, M);

for k = 1:M
    P1_R(:, :, k) = toRiemann(P1_samples(:, k).');
    P2_R(:, :, k) = toRiemann(P2_samples(:, k).');
    R1_nom(:, k) = toRiemann(P1_vals(k));
    R2_nom(:, k) = toRiemann(P2_vals(k));
end

%% Geodesic distance function
geodesic = @(r1, r2) 0.5 * real(acos(dot(r1, r2) ./ ...
                          (vecnorm(r1) .* vecnorm(r2))));

%% Initialize distances
wass_LP = zeros(1, M);
supp_dist = zeros(1, M);
lower_bound = zeros(1, M);

for k = 1:M
    fprintf("Solving LP at frequency point %d / %d...\n", k, M);

    % Compute cost matrix
    D = zeros(N, N);
    for i = 1:N
        for j = 1:N
            D(i, j) = geodesic(P1_R(:, i, k), P2_R(:, j, k));
        end
    end

    % Vectorize cost
    c = D(:);

    % Marginal constraints
    Aeq = zeros(2 * N, N^2);
    beq = ones(2 * N, 1) / N;

    % Row sums (P1)
    for i = 1:N
        Aeq(i, ((i-1)*N + 1):(i*N)) = 1;
    end

    % Column sums (P2)
    for j = 1:N
        Aeq(N + j, j:N:end) = 1;
    end

    % Solve LP
    f = c;
    options = optimoptions('linprog', 'Display', 'none');
    [gamma, fval] = linprog(f, [], [], Aeq, beq, zeros(N^2, 1), [], options);
    wass_LP(k) = fval;

    % Support distance
    supp_dist(k) = max(D(:));

    % Lower bound via nominal distance
    d_nom = geodesic(R1_nom(:, k), R2_nom(:, k));
    dev1 = mean(arrayfun(@(i) geodesic(P1_R(:, i, k), R1_nom(:, k)), 1:N));
    dev2 = mean(arrayfun(@(j) geodesic(P2_R(:, j, k), R2_nom(:, k)), 1:N));
    lower_bound(k) = abs(d_nom - dev1 - dev2);
end

%% Plotting
figure;
plot(omega, wass_LP, 'b-', 'LineWidth', 2); hold on;
plot(omega, supp_dist, 'r--', 'LineWidth', 2);
plot(omega, lower_bound, 'g-.', 'LineWidth', 2);
xlabel('Frequency \omega (rad/s)', 'FontSize', 12);
ylabel('Distance', 'FontSize', 12);
title('Frequency-Domain Distances via Linear Programming', 'FontSize', 14);
legend('Wasserstein-1 (LP)', 'Support Distance', 'Nominal Lower Bound');
grid on;