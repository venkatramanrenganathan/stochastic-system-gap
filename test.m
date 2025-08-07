% plot_two_systems_riemann_projection.m
% Shows 2 stochastic systems' supports, projections, and Riemann sphere embedding

clear; clc; close all;

%% Parameters
N = 100;  % Number of samples
spread = 0.05;

% Nominal points (distinct)
nom1 = 0.3 + 0.3j;
nom2 = -0.1 + 0.4j;

%% Generate complex samples for each system
samples1 = nom1 + spread * (randn(N,1) + 1j * randn(N,1));
samples2 = nom2 + spread * (randn(N,1) + 1j * randn(N,1));

x1 = real(samples1); y1 = imag(samples1);
x2 = real(samples2); y2 = imag(samples2);

%% Function: inverse stereographic projection
project_to_Riemann = @(z) [real(z) ./ (1 + abs(z).^2); ...
                           imag(z) ./ (1 + abs(z).^2); ...
                           abs(z).^2 ./ (1 + abs(z).^2)];

% Projections to Riemann sphere
proj1 = project_to_Riemann(samples1);
proj2 = project_to_Riemann(samples2);

%% Convex support sets (approximate as covariance ellipses)
ellipse_func = @(x, y) ...
    deal( ...
        mean([x y]), ...
        cov([x y]), ...
        @(mu,S) deal( ...
            mu(1) + cos(linspace(0,2*pi,200)) * sqrt(S(1,1)), ...
            mu(2) + sin(linspace(0,2*pi,200)) * sqrt(S(2,2)) ...
        ));

[mu1, Sigma1, ellipse1_func] = ellipse_func(x1, y1);
[mu2, Sigma2, ellipse2_func] = ellipse_func(x2, y2);
theta = linspace(0, 2*pi, 200);
[ell1_x, ell1_y] = ellipse1_func(mu1, Sigma1);
[ell2_x, ell2_y] = ellipse2_func(mu2, Sigma2);

% Complex ellipse points
ell1_cpx = ell1_x + 1j * ell1_y;
ell2_cpx = ell2_x + 1j * ell2_y;

% Project ellipse boundaries
proj1_ell = project_to_Riemann(ell1_cpx);
proj2_ell = project_to_Riemann(ell2_cpx);

% Project nominal points
proj_nom1 = project_to_Riemann(nom1);
proj_nom2 = project_to_Riemann(nom2);

%% Riemann sphere
[XS, YS, ZS] = sphere(50);
XS = XS * 0.5;
YS = YS * 0.5;
ZS = ZS * 0.5 + 0.5;

%% Plotting
figure('Color','w'); hold on; grid on; axis equal;
view(3); xlabel('Re'); ylabel('Im'); zlabel('Z');

% Complex plane support sets
fill3(ell1_x, ell1_y, zeros(size(ell1_x)), [0.9 0.9 1], ...
    'EdgeColor', 'b', 'LineWidth', 1.5, 'FaceAlpha', 0.5);
fill3(ell2_x, ell2_y, zeros(size(ell2_y)), [1 0.9 0.9], ...
    'EdgeColor', 'r', 'LineWidth', 1.5, 'FaceAlpha', 0.5);

% Samples
scatter3(x1, y1, zeros(N,1), 30, 'b', 'filled');
scatter3(x2, y2, zeros(N,1), 30, 'r', 'filled');

% Nominal points (on complex plane)
plot3(real(nom1), imag(nom1), 0, 'ko', 'MarkerFaceColor','b', 'MarkerSize', 6);
plot3(real(nom2), imag(nom2), 0, 'ko', 'MarkerFaceColor','r', 'MarkerSize', 6);

% Riemann sphere
surf(XS, YS, ZS, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.6 0.8 1]);

% Projection lines from samples
for i = 1:N
    plot3([x1(i), proj1(1,i)], [y1(i), proj1(2,i)], [0, proj1(3,i)], 'b:');
    plot3([x2(i), proj2(1,i)], [y2(i), proj2(2,i)], [0, proj2(3,i)], 'r:');
end

% Projection lines for nominal points
plot3([real(nom1), proj_nom1(1)], [imag(nom1), proj_nom1(2)], [0, proj_nom1(3)], 'b-', 'LineWidth', 2);
plot3([real(nom2), proj_nom2(1)], [imag(nom2), proj_nom2(2)], [0, proj_nom2(3)], 'r-', 'LineWidth', 2);

% Projected samples
scatter3(proj1(1,:), proj1(2,:), proj1(3,:), 30, 'b', 'filled');
scatter3(proj2(1,:), proj2(2,:), proj2(3,:), 30, 'r', 'filled');

% Projected support sets
plot3(proj1_ell(1,:), proj1_ell(2,:), proj1_ell(3,:), 'b-', 'LineWidth', 2);
plot3(proj2_ell(1,:), proj2_ell(2,:), proj2_ell(3,:), 'r-', 'LineWidth', 2);

% Projected nominal points
plot3(proj_nom1(1), proj_nom1(2), proj_nom1(3), 'ko', 'MarkerFaceColor', 'c', 'MarkerSize', 6);
plot3(proj_nom2(1), proj_nom2(2), proj_nom2(3), 'ko', 'MarkerFaceColor', 'm', 'MarkerSize', 6);

% Labels
text(real(nom1)+0.05, imag(nom1), 0.01, 'Nom 1', 'Color', 'b');
text(real(nom2)+0.05, imag(nom2), 0.01, 'Nom 2', 'Color', 'r');
plot3(0, 0, 0, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g'); % south pole
plot3(0, 0, 1, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g'); % north pole

title('Two Stochastic Systems in Complex Plane and Riemann Sphere');
legend({'Support 1 in \mathbb{C}', 'Support 2 in \mathbb{C}', ...
        'Samples 1', 'Samples 2', ...
        'Nominal 1', 'Nominal 2', 'Riemann sphere', ...
        'Projections 1', 'Projections 2', ...
        'Nom proj 1', 'Nom proj 2', ...
        'Projected samples 1', 'Projected samples 2', ...
        'Support 1 on sphere', 'Support 2 on sphere'}, ...
        'Location', 'northeastoutside');