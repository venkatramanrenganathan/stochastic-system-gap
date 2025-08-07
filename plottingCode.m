% plot_empirical_distribution_with_projection.m
% Empirical complex samples, inverse stereographic projections, Riemann sphere, and support

clear; clc; close all;

%% Parameters
N = 100;
nominal = 0.8 + 0.3j;
spread = 0.05;

%% Generate samples around nominal point
samples = nominal + spread * (randn(N, 1) + 1j * randn(N, 1));
x = real(samples);
y = imag(samples);

%% Compute ellipse (support set)
mu = [mean(x), mean(y)];
Sigma = cov([x y]);
theta = linspace(0, 2*pi, 200);
circle = [cos(theta); sin(theta)];
[U, S, ~] = svd(Sigma);
ellipse = U * sqrt(S) * circle;
ellipse(1,:) = ellipse(1,:) + mu(1);
ellipse(2,:) = ellipse(2,:) + mu(2);

%% Inverse stereographic projection
project_to_Riemann = @(z) [real(z) ./ (1 + abs(z).^2); ...
                           imag(z) ./ (1 + abs(z).^2); ...
                           abs(z).^2 ./ (1 + abs(z).^2)];

% Project samples to sphere (3 x N)
proj_R = project_to_Riemann(samples);
proj_x = proj_R(1:100,:);
proj_y = proj_R(101:200,:);
proj_z = proj_R(201:300,:);

% Project support ellipse
ellipse_cpx = ellipse(1,:) + 1j * ellipse(2,:);
ellipse_proj = project_to_Riemann(ellipse_cpx);
ellipse_proj_x = ellipse_proj(1,:);
ellipse_proj_y = ellipse_proj(2,:);
ellipse_proj_z = ellipse_proj(3,:);

% Project nominal point
nominal_proj = project_to_Riemann(nominal);

%% Riemann sphere
[XS, YS, ZS] = sphere(50); % unit diameter
XS = XS * 0.5;
YS = YS * 0.5;
ZS = ZS * 0.5 + 0.5;

%% Plot
figure('Color','w');
hold on; grid on; axis equal;
view(3);
xlabel('Re'); ylabel('Im'); zlabel('Z');

% Complex plane and support
fill3(ellipse(1,:), ellipse(2,:), zeros(size(ellipse(1,:))), ...
      [0.9 0.9 1], 'EdgeColor', 'b', 'LineWidth', 1.5, 'FaceAlpha', 0.5);
scatter3(x, y, zeros(N,1), 30, 'r', 'filled');

% Riemann sphere
surf(XS, YS, ZS, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.6 0.8 1]);

% Nominal and south pole
plot3(real(nominal), imag(nominal), 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
plot3(0, 0, 0, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
plot3(0, 0, 1, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');

% Projection lines from samples to sphere
for i = 1:N
    plot3([x(i), proj_x(i)], [y(i), proj_y(i)], [0, proj_z(i)], ...
        'k:', 'LineWidth', 1);
end

% Projection line for nominal point
plot3([real(nominal), nominal_proj(1)], ...
      [imag(nominal), nominal_proj(2)], ...
      [0, nominal_proj(3)], 'k-', 'LineWidth', 2);

% Projected sample points
scatter3(proj_x, proj_y, proj_z, 30, 'm', 'filled');

% Projected support set on sphere
plot3(ellipse_proj_x, ellipse_proj_y, ellipse_proj_z, 'm-', 'LineWidth', 2);

% Projected nominal point
plot3(nominal_proj(1), nominal_proj(2), nominal_proj(3), ...
      'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);

% Labels
text(real(nominal)+0.05, imag(nominal), 0.02, 'Nominal point', 'FontSize', 10);
text(0, 0, 1.05, 'North pole', 'HorizontalAlignment','center');
text(0, 0, -0.05, 'South pole', 'HorizontalAlignment','center');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 2);
set(a, 'linewidth', 2);
set(a, 'FontSize', 30);
set(gca,'fontweight','bold');
legend({'Support set in \mathbb{C}', 'Empirical samples', ...
        'Riemann sphere', 'Nominal (complex)', 'South pole', ...
        'Projection lines', 'Nominal projection', ...
        'Projected samples', 'Projected support'}, ...
        'Location', 'northeastoutside');