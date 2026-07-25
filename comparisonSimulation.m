%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code compares frequency-domain vs time-domain distances between two 
% stochastic systems.
% - Frequency domain distance is computed as Wasserstein-1 distance via 
% chordal metric on Riemann sphere followed by taking sup over omega
% - Time domain distance is computed as Wasserstein-1 distance via gap metric
%
% Copyrights @ 2026
% 
% Authors: Venkatraman Renganathan 
%          IIT Hyderabad, India
%
% Email: venkatraman@ai.iith.ac.in
%
% Date last updated: 25 July, 2026.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For reproducibility
rng(12);

fprintf('Generating data for two stochastic LTI systems \n');

% Define simulation parameters 
n = 2;                         % state dimension
m = 1;                         % control input dimension
l = 1;                         % output dimension
N = 100;                       % samples per system (modest: O(N^2) LPs)
p = 3;                         % dimension of parameter theta
w = logspace(-2, 2, 100);      % frequency grid (rad/s)
perturbationScale = 0.15;    % scales random directions 

% Theta follows Gaussian distribution parameters for both systems
meanTheta1 = [ 0.10; -0.05; 0.02];
meanTheta2 = [-0.08;  0.06; -0.01];
covarianceTheta1 = diag([0.15, 0.20, 0.10]).^2;        
covarianceTheta2 = diag([0.12, 0.18, 0.08]).^2;        

% Generate distinct nominal plant models for both the systems 
% System 1: second-order underdamped system
wn1 = 1.8; 
z1 = 0.35;
A1_0 = [0 1; -wn1^2 -2*z1*wn1]; 
B1_0 = [0; 1]; 
C1_0 = [1 0]; 
D1_0 = 0;

% System 2: second-order underdamped system - different resonance/damping
wn2 = 1.2; 
z2 = 0.55;
A2_0 = [0 1; -wn2^2 -2*z2*wn2]; 
B2_0 = [0; 1]; 
C2_0 = [1 0]; 
D2_0 = 0;

% Form the nominal state space models for both the systems
nominalModelSystem1 = ss(A1_0,B1_0,C1_0,D1_0);
nominalModelSystem2 = ss(A2_0,B2_0,C2_0,D2_0);

% Compute the gap between nominal models of both system 
nominalModelsGapDistance = gapmetric(nominalModelSystem1, nominalModelSystem2);

% Generate affine perturbation direction matrices for both systems
[A1_j, B1_j, C1_j] = generateRandomDirectionData(n,m,l,p,perturbationScale);
[A2_j, B2_j, C2_j] = generateRandomDirectionData(n,m,l,p,perturbationScale);

% Generate sample of theta parameters from Gaussian distribution 
sqrtCovarianceTheta1 = chol(covarianceTheta1, 'lower'); 
sqrtCovarianceTheta2 = chol(covarianceTheta2, 'lower'); 
%theta1 = meanTheta1 + sqrtCovarianceTheta1*randn(p,N);
%theta2 = meanTheta2 + sqrtCovarianceTheta2*randn(p,N);

theta1 = -0.1 + (0.1 - (-0.1)).*rand(p,N);
theta2 = -0.1 + (0.1 - (-0.1)).*rand(p,N);

% Build perturbed yet stable SS plant models for both systems 
perturbedModelsSystem1 = cell(1,N); 
perturbedModelsSystem2 = cell(1,N);
for k = 1:N
    [A1,B1,C1] = affinePerturbation(A1_0,B1_0,C1_0, A1_j,B1_j,C1_j, theta1(:,k));
    [A2,B2,C2] = affinePerturbation(A2_0,B2_0,C2_0, A2_j,B2_j,C2_j, theta2(:,k));
    % Form state space perturbed plant models using perturbed system data
    perturbedModelsSystem1{k} = ss(A1,B1,C1,D1_0);
    perturbedModelsSystem2{k} = ss(A2,B2,C2,D2_0);
end

%% ========== Compute the time domain distance ============================

fprintf('Computing Time Domain Distance Between Stochastic LTI Systems \n');

% Placeholder to store gap metric distance between models of both systems
gapMetricDistances = zeros(N,N);

% Get pairwise gap metric distance between perturbed models of both systems
for i = 1:N
    for j = 1:N
        gapMetricDistances(i,j) = gapmetric(perturbedModelsSystem1{i}, perturbedModelsSystem2{j});   
    end
end

% With gapMetricDistances as transport cost, find the time domain distance
timeDomainDistance = computeWassersteinLinearProgram(gapMetricDistances);

%% ========== Compute the frequency domain distance =======================

fprintf('Computing Frequency Domain Distance Between Stochastic LTI Systems \n');

% Placeholder to store the pointwise frequency distance at every frequency
frequencyWassersteinDistance = zeros(1, numel(w));

% At every frequency, compute the Wasserstein-1 distance
for iw = 1:numel(w)

    % Placeholders to store the frequency response data for both systems
    G1 = zeros(N,1); 
    G2 = zeros(N,1);
    
    % For every sample, get the frequency response data for both systems
    for k = 1:N
        G1(k) = squeeze(freqresp(perturbedModelsSystem1{k}, w(iw))); 
        G2(k) = squeeze(freqresp(perturbedModelsSystem2{k}, w(iw)));
    end

    % Map frequency responses to sphere for each systems
    R1 = projectToRiemannSphere(G1);        
    R2 = projectToRiemannSphere(G2);

    % Compute chordal distance between N x N points in Riemann sphere
    chordalDistances = pdist2(R1.', R2.', 'euclidean');   

    % Compute Wasserstein-1 distance with equal weights via LP
    frequencyWassersteinDistance(iw) = computeWassersteinLinearProgram(chordalDistances);
end

% Compute frequency domain distance as sup over frequency of frequencyWassersteinDistance
frequencyDomainDistance = max(frequencyWassersteinDistance);

%% ==================== Report the findings ==============================
fprintf('==== Reporting Distance Between Stochastic LTI Systems ====\n');
fprintf('Frequency-domain Distance  : %.6f\n', frequencyDomainDistance);
fprintf('Time-domain Distance       : %.6f\n', timeDomainDistance);
fprintf('Gap between nominal models : %.6f\n', nominalModelsGapDistance);
if(frequencyDomainDistance <= timeDomainDistance)
    fprintf('Frequency-domain distance <= time-domain distance. \n');
else
    fprintf('Frequency-domain distance > time-domain distance. \n');
end