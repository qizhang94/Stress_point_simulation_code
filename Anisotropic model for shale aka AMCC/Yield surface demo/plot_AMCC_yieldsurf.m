clear; clearvars -global; clc; close all;
addpath('C:\export_fig-master');

tiledlayout(1,3)
nexttile();
set(gcf, 'Color', 'w');

fimplicit(@(pminus, q) AMCC_simple(pminus, q, 1, 1), [0 10 0 11], 'b-', 'LineWidth', 1); hold on;
fimplicit(@(pminus, q) AMCC_simple(pminus, q, 1, 1.2), [0 10 0 11],'b:', 'LineWidth', 1);
fimplicit(@(pminus, q) AMCC_simple(pminus, q, 1, 0.8), [0 10 0 11], 'b-.', 'LineWidth', 1);
axis equal;

xlabel('$-p$', 'Interpreter','latex');
ylabel('$q$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',16); % Adjust label size
set(get(gca,'YLabel'),'FontSize',16);
lgd = legend('$\beta_p = 1$','$\beta_p = 1.2$','$\beta_p = 0.8$','location', 'best');
legend boxoff
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
% lgd.Layout.Tile = 0; % assign tile location


nexttile();
set(gcf, 'Color', 'w');

fimplicit(@(pminus, q) AMCC_simple(pminus, q, 1, 1), [0 10 0 11], 'r-', 'LineWidth', 1); hold on;
fimplicit(@(pminus, q) AMCC_simple(pminus, q, 1.2, 1), [0 10 0 11],'r:', 'LineWidth', 1);
fimplicit(@(pminus, q) AMCC_simple(pminus, q, 0.9, 1), [0 10 0 11], 'r-.', 'LineWidth', 1);
axis equal;

xlabel('$-p$', 'Interpreter','latex');
ylabel('$q$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',16); % Adjust label size
set(get(gca,'YLabel'),'FontSize',16);
lgd = legend('$\alpha_p = 1$','$\alpha_p = 1.2$','$\alpha_p = 0.9$','location', 'best');
legend boxoff
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
% lgd.Layout.Tile = 2; % assign tile location

nexttile();
set(gcf, 'Color', 'w');
fimplicit(@(pminus, q) AMCC_simple(pminus, q, 1, 1), [0 10 0 11], 'k-', 'LineWidth', 1); hold on;
axis equal;

xlabel('$-p^*$', 'Interpreter','latex');
ylabel('$q^*$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',16); % Adjust label size
set(get(gca,'YLabel'),'FontSize',16);


%export_fig AMCC_pq_surface.png -m7 -painters;

function y = AMCC_simple(pminus, q, alpha_p, beta_p)

p = -pminus; % p = (sigma_1 + 2*sigma_3)/3; q > 0 = sigma_3 - sigma_1 (both sigma_3 and sigma_1 are negative)
n = [0; 0; 1];
m = n*n';
c1p = 0.7; c2p = -0.36; c3p = 0.6; % 1, 0, 0 Isotropic
% alpha_p = c1p + c2p + c3p;
% beta_p = c1p;
% gamma_p = c1p + c3p/2; % fix at 1

sigma1 = p - 2/3*q; % in z direction
sigma3 = p + 1/3*q; % in x and y
sigma2 = sigma3;

sigma1_fic = sigma1*alpha_p;
sigma2_fic = sigma2*beta_p;
sigma3_fic = sigma3*beta_p;

p_fic = (sigma1_fic + sigma2_fic + sigma3_fic)/3;
q_fic = abs(sigma1_fic - sigma3_fic);

M = 1; % CSL slope
pc = -10; % MPa

y = q_fic.^2/M^2 + p.*(p-pc);

end


