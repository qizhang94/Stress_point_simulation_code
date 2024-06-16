clear; clearvars -global; clc; close all;
restoredefaultpath;
export = false;
format short
% if export
%     addpath('C:\export_fig-master');
% end

% Investigate the trend of triaxial undrained test
K = 17.5*1e3; % kPa
G = 10.5*1e3; % kPa
E = 9*K*G/(3*K + G);
nu = (3*K - 2*G)/2/(3*K + G);
Props = [E, nu, 0, 0, 10.01]; % unit of c is kPa

pconf0 = 400; % Initial confining pressure (isotropic) [kPa]
stress_old = [-pconf0; -pconf0; -pconf0; 0; 0; 0];
stress_new = stress_old;

hsv_old = 0; % History variables: PEEQ
hsv_new = hsv_old;
strain_accu = zeros(6,1); % Strain change "from" the isotropic compression state

comp_incr = -1e-3; % Compression strain amount increment
ntolstep = 100; % Number of total simulation steps

q_data = [0]; ea_data = [0]; p_data = [pconf0]; % Resulting data for plot

for i = 1:1:ntolstep
    strain_iter = [-comp_incr/2; -comp_incr/2; comp_incr; 0; 0; 0]; % no volumetric strain
   [stress_new,hsv_new,DDSDDE]=DP_UMAT(Props,stress_old,strain_iter,hsv_old);

    % Save state
    stress_old = stress_new;
    hsv_old = hsv_new;
    strain_accu = strain_accu + strain_iter;

    q_data = [q_data; stress_old(1) - stress_old(3)];
    ea_data = [ea_data; -strain_accu(3)];
    p_data = [p_data; -sum(stress_old(1:3))/3];

end

figure(1);
set(gcf, 'Color', 'w');
fig1 = plot(ea_data,  q_data, 'b', 'LineWidth', 1); grid on;
xlabel('Compressive axial strain $\varepsilon_a$', 'Interpreter','latex'); 
ylabel('von Mises stress $q = \|\sigma_1 - \sigma_3\|$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
% if export
%     export_fig(sprintf('../Cropped_PDF_figures/CUeaq_OCR%d.png', -pc0/pconf0), '-m7', '-painters');
% end

figure(2);
set(gcf, 'Color', 'w');
PWP = q_data/3+ pconf0 - p_data;
plot(ea_data, PWP, 'b', 'LineWidth', 1); grid on;
set(gca, 'YDir','reverse');
xlabel('Compressive axial strain $\varepsilon_a$', 'Interpreter','latex'); 
ylabel('PWP (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
% if export
%     export_fig(sprintf('../Cropped_PDF_figures/CUeapwp_OCR%d.png', -pc0/pconf0), '-m7', '-painters');
% end

% figure(3);
% set(gcf, 'Color', 'w');
% semilogx(p_data, v_data, 'g', 'LineWidth', 1); grid on; hold on;
% semilogx([p_data; -pc0], Props(4) - Props(2)*log([p_data; -pc0]), 'b-.', 'LineWidth', 1);
% semilogx([p_data; -pc0], Props(4)  - (Props(2)-Props(3))*log(2) - Props(2)*log([p_data; -pc0]), 'r-.', 'LineWidth', 1);
% semilogx([pconf0, -pc0],[v_start, Props(4) - Props(2)*log(-pc0)], 'k--');
% lgd = legend('Simulated CU', 'ICL, $\lambda$ line', 'CSL', 'Swelling $\kappa$ line');
% lgd.FontSize = 14;
% lgd.Interpreter = 'Latex';
% xlabel('$-p''$ (kPa)', 'Interpreter','latex');
% ylabel('Specific volume $\bar{v}$', 'Interpreter','latex');
% set(gca, 'TickLabelInterpreter','Latex');
% set(gca, 'FontSize', 13); % Adjust axis tick number/value size
% set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
% set(get(gca,'YLabel'),'FontSize',15);
% if export
%     export_fig(sprintf('../Cropped_PDF_figures/CUvlogp_OCR%d.png', -pc0/pconf0), '-m7', '-painters');
% end

figure(4);
set(gcf, 'Color', 'w');
plot(-p_data, q_data, 'g', 'LineWidth', 1); grid on; hold on; axis equal; % Slope = 3
xlabel('$p''$ (kPa)', 'Interpreter','latex');
ylabel('von Mises stress $q = \|\sigma_1 - \sigma_3\|$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
% if export
%     export_fig(sprintf('../Cropped_PDF_figures/CUpq_OCR%d.png', -pc0/pconf0), '-m7', '-painters');
% end
% 
% function [] = MCC_ys(pc, M)  % Add MCC yield surface plot
%     fplot(@(x) sqrt(x.*(-pc - x))*M, [0, -pc], 'k--'); hold on
% end

