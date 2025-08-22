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

% Set parameters for frequency domain simulation 
numSamples = 100;       % # of frequency response samples at every frequency
numFrequencies = 1000;  % # of frequency points

% Create discretised frequency space between [0.1, 100] with M points
omegaMin = 0.1;                                       % Minimum frequency
omegaMax = 1000;                                      % Maximum frequency 
Omega = linspace(omegaMin, omegaMax, numFrequencies); % discretised frequency space

% Define nominal plant models for two different systems P1 & P2
P1NominalModel = @(w) 1 ./ (1 + 0.5j * w);
P2NominalModel = @(w) 1 ./ ((1 + 0.2j * w) .* (1 + 0.7j * w));

% Evaluate the nominal models of both the systems at all frequencies 
P1NominalResponses = P1NominalModel(Omega);
P2NominalResponses = P2NominalModel(Omega);

% Placeholders to store N frequency response samples at each of the M frequency points
P1Samples = zeros(numSamples, numFrequencies);
P2Samples = zeros(numSamples, numFrequencies);

% Define the perturbation scale
perturbationFactor = 0.01; 

% Generate frequency response samples
for i = 1:numSamples
    % Get a random complex perturbation for every frequency for both systems
    ithSampleNoiseSystem1 = perturbationFactor * (randn(1, numFrequencies) + 1j * randn(1, numFrequencies));
    ithSampleNoiseSystem2 = perturbationFactor * (randn(1, numFrequencies) + 1j * randn(1, numFrequencies));
    % Generate the ith sample at every frequency by combining the nominal  
    % response & the random perturbation
    P1Samples(i, :) = P1NominalResponses .* (1 + ithSampleNoiseSystem1);
    P2Samples(i, :) = P2NominalResponses .* (1 + ithSampleNoiseSystem2);
end

% Placeholders for storing inverse stereographic projected quantities
empiricalDistributionP1R = zeros(3, numSamples, numFrequencies);
empiricalDistributionP2R = zeros(3, numSamples, numFrequencies);
R1NominalValues = zeros(3, numFrequencies);
R2NominalValues = zeros(3, numFrequencies);

% Projection items from complex plane to Riemann sphere (R) for both systems
for k = 1:numFrequencies
    % Project distribution of frequency responses in complex field to R
    empiricalDistributionP1R(:, :, k) = inverseStereographicProjection(P1Samples(:, k).');
    empiricalDistributionP2R(:, :, k) = inverseStereographicProjection(P2Samples(:, k).');
    % Project nominal frequency responses from complex field to R
    R1NominalValues(:, k) = inverseStereographicProjection(P1NominalResponses(k));
    R2NominalValues(:, k) = inverseStereographicProjection(P2NominalResponses(k));
end

% Place holders to store the computed distance quantities
type1WassersteinDistance = zeros(1, numFrequencies);
upperBoundViaSupportDistance = zeros(1, numFrequencies);
triangleInequalityDistanceLowerBound = zeros(1, numFrequencies);

% Loop for every frequency & compute the distance & its lower & upper bounds
for k = 1:numFrequencies
    
    % Display that we are solving the linear program at frequency k
    fprintf("Solving LP at frequency %d / %d...\n", k, numFrequencies);

    % Compute travel cost matrix
    D = zeros(numSamples, numSamples);
    for i = 1:numSamples
        for j = 1:numSamples
            D(i, j) = computeDistance(empiricalDistributionP1R(:, i, k), empiricalDistributionP2R(:, j, k));
        end
    end

    % Vectorize the cost
    cost = D(:);

    % Marginal constraints
    Aeq = zeros(2 * numSamples, numSamples^2);
    beq = ones(2 * numSamples, 1) / numSamples;

    % Row sums (P1)
    for i = 1:numSamples
        Aeq(i, ((i-1)*numSamples + 1):(i*numSamples)) = 1;
    end

    % Column sums (P2)
    for j = 1:numSamples
        Aeq(numSamples + j, j:numSamples:end) = 1;
    end

    % Solve the Linear Program
    options = optimoptions('linprog', 'Display', 'none');
    [gamma, fval] = linprog(cost, [], [], Aeq, beq, zeros(numSamples^2, 1), [], options);

    % Get the Wasserstein distance from the solution of LP
    type1WassersteinDistance(k) = fval;

    % % Directed distance
    % hAB = max(min(D, [], 2));  % max over A of min over B
    % hBA = max(min(D, [], 1));  % max over B of min over A

    % Support Hausdorff distance
    upperBoundViaSupportDistance(k) = max(D(:));
    %upperBoundViaSupportDistance(k) = max(hAB, hBA);

    % Compute the distance between nominal models at frequency k
    nominalModelsDistance = computeDistance(R1NominalValues(:, k), R2NominalValues(:, k));

    % Deviation of first system from its nominal frequency response
    system1DeviationFromNominal = mean(arrayfun(@(i) computeDistance(empiricalDistributionP1R(:, i, k), R1NominalValues(:, k)), 1:numSamples));

    % Deviation of second system from its nominal frequency response
    system2DeviationFromNominal = mean(arrayfun(@(j) computeDistance(empiricalDistributionP2R(:, j, k), R2NominalValues(:, k)), 1:numSamples));

    % Compute distance lower bound via triangle inequality due to deviation
    % with q = 1
    triangleInequalityDistanceLowerBound(k) = max(nominalModelsDistance - system1DeviationFromNominal - system2DeviationFromNominal, 0);

end

%% Display Results
fprintf('Distance Upper Bound = %.4f\n', max(upperBoundViaSupportDistance));
fprintf('Type-1 Wasserstein Distance = %.4f\n', max(type1WassersteinDistance));
fprintf('Distance Lower Bound = %.4f\n', max(triangleInequalityDistanceLowerBound));

%% Plotting Code

% Plot the type-q Gap-Wasssertein distance, its lower & upper bounds
figure;
slx1 = semilogx(Omega, type1WassersteinDistance, 'b-.'); hold on;
slx2 = semilogx(Omega, upperBoundViaSupportDistance, 'r-.');
slx3 = semilogx(Omega, triangleInequalityDistanceLowerBound, '-.');
slx3.Color = 'magenta';
yl1 = yline(max(type1WassersteinDistance),'b-','$\mathrm{d}_{1}(P_1, P_2)$', 'Interpreter','latex');
yl2 = yline(max(upperBoundViaSupportDistance),'r-','$\mathrm{d^{\mathrm{R}}_{sup}}(P_1, P_2)$', 'Interpreter','latex');
yl3 = yline(max(triangleInequalityDistanceLowerBound),'m-','Greatest Lower Bound', 'Interpreter','latex');
yl1.LabelHorizontalAlignment = 'center';
yl2.LabelHorizontalAlignment = 'left';
yl3.LabelHorizontalAlignment = 'left';
yl1.FontSize = 30;
yl2.FontSize = 30;
yl3.FontSize = 30;
yl1.LineWidth = 8;
yl2.LineWidth = 8;
yl3.LineWidth = 8;
xlabel('frequency $\omega$ (rad/s)');
ylabel('distance');
legend('$W^{1}_{1}\left(P_{R_{1}}(\omega), P_{R_{2}}(\omega)\right)$', '$\mathrm{d^{\mathrm{R}}_{sup}}(P_1, P_2, \omega)$', 'Lower Bound at $\omega$', 'Location','southeast');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 8);
set(a, 'linewidth', 8);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');