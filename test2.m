function stochastic_system_distance_demo
% Compares frequency-domain vs time-domain distances between two stochastic systems.
% - Frequency domain: Wasserstein-1 using CHORDAL (Euclidean) metric on Riemann sphere; sup over omega
% - Time domain: Wasserstein-1 using ν-gap metric (gapmetric)
%
% Outputs printed at the end:
%   W1_freq_sup  : sup_omega W1_chordal(omega)
%   W1_time_gap  : W1 with ν-gap as transport cost (sample clouds of plants)
%   nugap_nom    : ν-gap between nominal plants (sanity check)
%
% Requires Optimization Toolbox (linprog). Robust Control Toolbox (gapmetric) is optional.


clear all; close all; clc;

rng(12);                                % reproducibility

%% ------------------ User parameters ------------------
N = 100;                                  % samples per system (keep modest: O(N^2) LPs)
p = 3;                                   % parameter dimension
w = logspace(-2, 2, 100);                 % frequency grid (rad/s)
perturb_scale = 0.15;                    % scales random directions; keep small to ensure stability

% Parameter Gaussian specs
m1 = [ 0.10; -0.05; 0.02];
m2 = [-0.08;  0.06; -0.01];
S1 = diag([0.15, 0.2, 0.1]).^2;          % covariance for theta1
S2 = diag([0.12, 0.18, 0.08]).^2;        % covariance for theta2

%% ------------------ Nominal (distinct) SISO systems ------------------
% System 1: second-order underdamped
wn1 = 1.8; z1 = 0.35;
A1_0 = [0 1; -wn1^2 -2*z1*wn1]; B1_0 = [0; 1]; C1_0 = [1 0]; D1_0 = 0;

% System 2: different resonance/damping
wn2 = 1.2; z2 = 0.55;
A2_0 = [0 1; -wn2^2 -2*z2*wn2]; B2_0 = [0; 1]; C2_0 = [1 0]; D2_0 = 0;

sys1_nom = ss(A1_0,B1_0,C1_0,D1_0);
sys2_nom = ss(A2_0,B2_0,C2_0,D2_0);

%% ------------------ Affine perturbation directions ------------------
[n, m, l] = deal(2,1,1);
[A1_j, B1_j, C1_j] = rand_dirs(n,m,l,p,perturb_scale);
[A2_j, B2_j, C2_j] = rand_dirs(n,m,l,p,perturb_scale);

%% ------------------ Sample parameters (Gaussian) ------------------
L1 = chol(S1, 'lower'); theta1 = m1 + L1*randn(p,N);
L2 = chol(S2, 'lower'); theta2 = m2 + L2*randn(p,N);

%% ------------------ Build perturbed SS models (ensure stability) ------------------
sys1_set = cell(1,N); sys2_set = cell(1,N);
for k = 1:N
    [A1,B1,C1] = affine_perturb(A1_0,B1_0,C1_0, A1_j,B1_j,C1_j, theta1(:,k));
    [A2,B2,C2] = affine_perturb(A2_0,B2_0,C2_0, A2_j,B2_j,C2_j, theta2(:,k));
    A1 = ensure_hurwitz(A1); A2 = ensure_hurwitz(A2);
    sys1_set{k} = ss(A1,B1,C1,D1_0);
    sys2_set{k} = ss(A2,B2,C2,D2_0);
end

%% ------------------ FREQUENCY-DOMAIN: W1_chordal (sup over omega) ------------------
% For each omega, compute stochastic frequency responses; chordal distance via Riemann sphere (R=1/2)
W1_freq = zeros(1, numel(w));

for iw = 1:numel(w)
    % FRF samples
    G1 = zeros(N,1); G2 = zeros(N,1);
    for k = 1:N
        G1(k) = squeeze(freqresp(sys1_set{k}, w(iw))); % complex scalar
        G2(k) = squeeze(freqresp(sys2_set{k}, w(iw)));
    end

    % Map to sphere and compute chordal (Euclidean in 3D) distances
    R1 = toSphere(G1);        % 3 x N points on sphere of radius 1/2 centered at (0,0,1/2)
    R2 = toSphere(G2);

    D = pdist2(R1.', R2.', 'euclidean');   % chordal distance matrix N x N

    % Wasserstein-1 with equal weights via LP
    W1_freq(iw) = wasserstein_lp(D);
end

W1_freq_sup = max(W1_freq);

%% ------------------ TIME-DOMAIN: W1 with ν-gap as cost ------------------

Dgap = zeros(N,N);
for i = 1:N
    for j = 1:N
        Dgap(i,j) = gapmetric(sys1_set{i}, sys2_set{j});   % ν-gap in [0,1]
    end
end

W1_time_gap = wasserstein_lp(Dgap);

%% ------------------ (Optional) ν-gap between nominals ------------------
nugap_nom = gapmetric(sys1_nom, sys2_nom);

%% ------------------ REPORT ------------------
fprintf('\n===== Stochastic System Distance Report =====\n');
fprintf('Samples per system    : N = %d\n', N);
fprintf('Frequency grid points : M = %d\n', numel(w));
fprintf('Parameter dims (p)    : %d\n', p);
fprintf('--- Results ---\n');
fprintf('Frequency-domain  W1 (chordal), sup over omega : %.6f\n', W1_freq_sup);
fprintf('Time-domain       W1 (ν-gap cost)               : %.6f\n', W1_time_gap);
fprintf('ν-gap between nominal plants                    : %.6f\n', nugap_nom);
if(W1_freq_sup > W1_time_gap + 1e-8)
    fprintf('NOTE: The frequency-domain value exceeds the time-domain ν-gap Wasserstein.\n');
    fprintf('      This can happen; ν-gap ≤ (graph gap), while chordal-W1 is not ν-gap.\n');
end
fprintf('=============================================\n');

end % main

%% ---------- Helpers ----------

function [A,B,C] = affine_perturb(A0,B0,C0, Aj,Bj,Cj, theta)
    A = A0; B = B0; C = C0;
    for r = 1:numel(theta)
        A = A + theta(r)*Aj(:,:,r);
        B = B + theta(r)*Bj(:,:,r);
        C = C + theta(r)*Cj(:,:,r);
    end
end

function A = ensure_hurwitz(A)
    % Shift any right-half-plane eigenvalues to ensure stability.
    ev = eig(A);
    maxRe = max(real(ev));
    if maxRe >= -0.2
        A = A - (maxRe + 0.3)*eye(size(A));
    end
end

function [Aj,Bj,Cj] = rand_dirs(n,m,l,p,scale)
    Aj = zeros(n,n,p);
    Bj = zeros(n,m,p);
    Cj = zeros(l,n,p);
    for r = 1:p
        Aj(:,:,r) = scale * (randn(n)*0.5);
        Bj(:,:,r) = scale * (randn(n,m)*0.5);
        Cj(:,:,r) = scale * (randn(l,n)*0.5);
    end
end

function R = toSphere(z)
    % Inverse stereographic projection onto sphere radius 1/2, center (0,0,1/2)
    % z: N x 1 complex; returns 3 x N matrix of 3D points
    x = real(z(:)).'; y = imag(z(:)).'; r2 = x.^2 + y.^2;
    X = x ./ (1 + r2);
    Y = y ./ (1 + r2);
    Z = r2 ./ (1 + r2);
    R = [X; Y; Z]; % points on sphere of radius 1/2 centered at (0,0,1/2)
end

function val = wasserstein_lp(D)
    % Solve W1 with uniform marginals via LP: min <D, Gamma>
    % subject to Gamma*1 = 1/N, Gamma^T*1 = 1/N, Gamma >= 0
    [N1,N2] = size(D);
    assert(N1==N2, 'Use equal sample counts on both sides.');
    N = N1;
    f = D(:);

    % Equality constraints
    Aeq = zeros(2*N, N*N);
    beq = ones(2*N,1) / N;

    % Row sums (source)
    for i = 1:N
        idx = (i-1)*N + (1:N);
        Aeq(i, idx) = 1;
    end
    % Column sums (target)
    for j = 1:N
        Aeq(N+j, j:N:end) = 1;
    end

    lb = zeros(N*N,1);
    opts = optimoptions('linprog','Display','none','Algorithm','dual-simplex');
    [gamma, fval, exitflag] = linprog(f, [], [], Aeq, beq, lb, [], opts);
    if exitflag <= 0
        error('linprog failed (exitflag=%d). Consider reducing N.', exitflag);
    end
    val = fval; % Wasserstein-1 cost with given metric
end