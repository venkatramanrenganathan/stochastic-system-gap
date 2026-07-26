%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency-domain stochastic distance
%
% Type-2 Wasserstein-induced distance on the Riemann sphere
%
% Copyrights @ 2026
% 
% Authors: Venkatraman Renganathan 
%          IIT Hyderabad, India
%
% Email: venkatraman@ai.iith.ac.in
%
% Date last updated: 26 July, 2026.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

%% Plotting settings

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Set parameters for frequency domain simulation 
numSamples = 100;       % # of frequency response samples at every frequency
numFrequencies = 1000;  % # of frequency points

% Create discretised frequency space between [0.1, 1000] with numFrequencies points
Omega = linspace(0.1,1000,numFrequencies);

% Define nominal plant models for two different systems P1 & P2
P1NominalModel = @(w) 1 ./ (1 + 0.5j*w);
P2NominalModel = @(w) 1 ./ ((1 + 0.2j*w).*(1 + 0.7j*w));

% Evaluate the nominal models of both the systems at all frequencies 
P1NominalResponses = P1NominalModel(Omega);
P2NominalResponses = P2NominalModel(Omega);

% Placeholders to store N frequency response samples at each of the M frequency points
P1Samples = zeros(numSamples,numFrequencies);
P2Samples = zeros(numSamples,numFrequencies);

% Define the perturbation scale
perturbationFactor = 0.01;

% Generate frequency response samples
for i = 1:numSamples
     % Get a random complex perturbation for every frequency for both systems
    ithSampleNoiseSystem1 = perturbationFactor * (rand(1,numFrequencies) + 1j*rand(1,numFrequencies));
    ithSampleNoiseSystem2 = perturbationFactor * (rand(1,numFrequencies) + 1j*rand(1,numFrequencies));
    % Generate the ith sample at every frequency by combining the nominal  
    % response & the random perturbation
    P1Samples(i,:) = P1NominalResponses .* (1 + ithSampleNoiseSystem1);
    P2Samples(i,:) = P2NominalResponses .* (1 + ithSampleNoiseSystem2);
end

% Placeholders for storing inverse stereographic projected quantities
empiricalDistributionP1R = zeros(3, numSamples, numFrequencies);
empiricalDistributionP2R = zeros(3, numSamples, numFrequencies);
R1NominalValues = zeros(3, numFrequencies);
R2NominalValues = zeros(3, numFrequencies);

% Projection items from complex plane to Riemann sphere (R) for both systems
for k = 1:numFrequencies
    % Project distribution of frequency responses in complex field to R
    empiricalDistributionP1R(:,:,k) = inverseStereographicProjection(P1Samples(:,k).');
    empiricalDistributionP2R(:,:,k) = inverseStereographicProjection(P2Samples(:,k).');
    % Project nominal frequency responses from complex field to R
    R1NominalValues(:, k) = inverseStereographicProjection(P1NominalResponses(k));
    R2NominalValues(:, k) = inverseStereographicProjection(P2NominalResponses(k));
end

% Place holders to store the computed distance quantities
type2WassersteinDistance = zeros(1,numFrequencies);
upperBoundViaSupportDistance = zeros(1,numFrequencies);
triangleInequalityDistanceLowerBound = zeros(1,numFrequencies);

%% Linear programming options
options = optimoptions('linprog','Display','none');

% Equality constraints for transport plan
% Every empirical distribution has equal weights 1/N.
Aeq = zeros(2*numSamples,numSamples^2);
beq = ones(2*numSamples,1)/numSamples;

% Row constraints
for i = 1:numSamples
    Aeq(i, ((i-1)*numSamples + 1):(i*numSamples)) = 1;
end

% Column constraints
for k = 1:numSamples
    Aeq(numSamples+k, k:numSamples:end) = 1;
end

% Loop for every frequency & compute the distance & its lower & upper bounds
for k = 1:numFrequencies

    % Display that we are solving for frequency k
    fprintf("Solving for frequency %d / %d...\n", k, numFrequencies);
    
    % Compute travel cost matrix
    D = zeros(numSamples,numSamples);
    for i = 1:numSamples
        for j = 1:numSamples
            D(i,j) = computeDistance(empiricalDistributionP1R(:,i,k), empiricalDistributionP2R(:,j,k));
        end
    end

    % Type-2 Wasserstein distance. No square root taken as paper defines W_q^q.
    costType2 = D(:).^2;

    % Solve the Linear Program
    [~,fval,exitflag] = linprog(costType2, [],[], Aeq,beq, zeros(numSamples^2,1), [], options);

    if exitflag <= 0
        warning('Type-2 LP did not converge at frequency %d.',k);
    end

    % This is W_2^2 in conventional OT terminology & it is the q=2 distance used in the paper.
    type2WassersteinDistance(k) = fval;

    % get upper bound
    upperBoundViaSupportDistance(k) = max(D(:))^2;

    %% ---------------------------------------------------------------
    % Step 4: Nominal chordal distance
    % ---------------------------------------------------------------

    nominalModelsDistance = computeDistance(R1NominalValues(:,k), R2NominalValues(:,k));
    deviationSystem1 = zeros(numSamples,1);
    deviationSystem2 = zeros(numSamples,1);
    for i = 1:numSamples
        deviationSystem1(i) = computeDistance(empiricalDistributionP1R(:,i,k), R1NominalValues(:,k));
        deviationSystem2(i) = computeDistance(empiricalDistributionP2R(:,i,k), R2NominalValues(:,k));
    end

    % Matrix Delta_nom^(ij)
    DeltaNom = deviationSystem1 + deviationSystem2.';

    % Note that Delta_dev = d_chord(R1,R2) - sum_ij DeltaNom_ij pi_ij and
    % lower bound = inf_pi (Delta_dev)_+^2.
    % Since (x)_+^2 is monotonically increasing in x, minimizing (Delta_dev)_+^2 is equivalent to maximizing sum DeltaNom_ij pi_ij.
    % linprog minimizes, so maximize DeltaNom by minimizing -DeltaNom.
    costLowerBound = -DeltaNom(:);

    % Solve LP
    [~,fvalLower,exitflagLower] = linprog(costLowerBound, [],[], Aeq,beq, zeros(numSamples^2,1), [], options);

    % Maximum possible expected nominal deviation
    maximumExpectedDeviation = -fvalLower;

    % Delta_dev corresponding to the minimizing lower-bound expression
    DeltaDev = nominalModelsDistance - maximumExpectedDeviation;

    % Square as q = 2
    triangleInequalityDistanceLowerBound(k) = max(DeltaDev,0)^2;

end

%% ========================================================================
% Infer Global distance and bounds
% ========================================================================
frequencyDomainDistance = max(type2WassersteinDistance);
upperBound = max(upperBoundViaSupportDistance);
lowerBound = max(triangleInequalityDistanceLowerBound);

%% Display
fprintf('\n Type-2 distance = %.6f\n', frequencyDomainDistance);
fprintf('Upper bound    = %.6f\n', upperBound);
fprintf('Lower bound    = %.6f\n', lowerBound);

%% Plot Results
figure;
semilogx(Omega, type2WassersteinDistance, 'LineWidth',4);
hold on;
semilogx(Omega, upperBoundViaSupportDistance, 'LineWidth',4);
semilogx(Omega, triangleInequalityDistanceLowerBound, 'LineWidth',4);
yline(frequencyDomainDistance, '--b', 'LineWidth',3);
yline(upperBound, '--r', 'LineWidth',3);
yline(lowerBound, '--o', 'LineWidth',3);
xlabel('$\omega$ (rad/s)', 'Interpreter','latex');
ylabel('Distance', 'Interpreter','latex');
legend('$\widehat d_2(P_1,P_2,\omega)$', 'Upper bound', 'Lower bound', '$\widehat d_2(P_1,P_2)$', 'Upper bound', 'Lower bound', 'Location','best', 'Interpreter','latex');
grid on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 45);
set(gca,'fontweight','bold');