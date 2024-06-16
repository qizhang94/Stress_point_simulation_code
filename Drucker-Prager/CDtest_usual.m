clear; clearvars -global; clc; close all;
restoredefaultpath;
export = false;
format short

% Investigate the trend of triaxial drained test
K = 17.5*1e3; % kPa
G = 10.5*1e3; % kPa
E = 9*K*G/(3*K + G);
nu = (3*K - 2*G)/2/(3*K + G);
Props = [E, nu, 29*pi/180, 29*pi/180*1/3, 10.01]; % unit of c is kPa

pconf0 = 600; % Initial confining pressure (isotropic) [kPa]
stress_old = [-pconf0; -pconf0; -pconf0; 0; 0; 0];
stress_new = stress_old;

hsv_old = 0; % History variables: PEEQ
hsv_new = hsv_old;
strain_accu = zeros(6,1); % Strain change "from" the isotropic compression state

comp_incr = -2e-3; % Compression strain amount increment
ntolstep = 100; % Number of total simulation steps

b = [-pconf0; -pconf0; comp_incr; 0; 0; 0]; % RHS of the non-linear equation
q_data = [0]; ea_data = [0]; ev_data = [0]; p_data = [pconf0]; % Resulting data for plot

for i = 1:1:ntolstep
    strain_iter = [-nu*comp_incr; -nu*comp_incr; comp_incr; 0; 0; 0]; % assuming elastic deformation
    error_tol = 1e-12;
    k = 0; % Outer iteration counter, inner is in the UMAT function
    while 1 > 0
        [stress_new,hsv_new,DDSDDE]=DP_UMAT(Props,stress_old,strain_iter,hsv_old);
%         [stress_new,~,DDSDDE]=MohrCoulomb_UMAT_3D(Props,stress_old,strain_iter);
        R = stress_new - b; % Residual
        R(3) = 0; % Because we know the epsilon_zz, then the increment should be 0
        if k == 0
            R_ini = R + 1; % Initial residual for this step, the addition of 1 is arbitrary
        end
        % Output convergence hisotry
        disp(['load step = ', num2str(i),'; iter = ', num2str(k),...
            '; conf = ', num2str(stress_new(1)), '; axial = ', num2str(stress_new(3))]);

        if norm(R)/norm(R_ini) < error_tol || norm(R) < error_tol
            break;
        end
        if k > 50  % Not converge
            break;
        end

        DDSDDE(3,:) = 0; DDSDDE(:, 3) = 0; DDSDDE(3,3) = 1;
        strain_iter = strain_iter - DDSDDE\R; % Guarantee the third component is unchanged
        k = k + 1;
    end
    % Save state
    stress_old = stress_new;
    hsv_old = hsv_new;
    strain_accu = strain_accu + strain_iter;

    q_data = [q_data; stress_old(1) - stress_old(3)];
    ea_data = [ea_data; -strain_accu(3)];
    ev_data = [ev_data; -sum(strain_accu(1:3))];
    p_data = [p_data; -sum(stress_old(1:3))/3];

end

figure(1);
set(gcf, 'Color', 'w');
plot(ea_data,  q_data/1e3, 'b', 'LineWidth', 1); grid on;
xlabel('Compressive axial strain $\varepsilon_a$', 'Interpreter','latex'); 
ylabel('von Mises stress $q = \|\sigma_1 - \sigma_3\|$ (MPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


figure(2);
set(gcf, 'Color', 'w');
plot(ea_data, -ev_data, 'r', 'LineWidth', 1); grid on; axis equal
% set(gca, 'YDir','reverse');
xlabel('Compressive axial strain $\varepsilon_a$', 'Interpreter','latex');
ylabel('Volumetric strain $\varepsilon_v$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


% 
% 
% 
% figure(3);
% set(gcf, 'Color', 'w');
% semilogx(p_data, v_data, 'g', 'LineWidth', 1); grid on; hold on;
% semilogx([p_data; -pc0], Props(4) - Props(2)*log([p_data; -pc0]), 'b-.', 'LineWidth', 1);
% semilogx([p_data; -pc0], Props(4)  - (Props(2)-Props(3))*log(2) - Props(2)*log([p_data; -pc0]), 'r-.', 'LineWidth', 1);
% semilogx([pconf0, -pc0],[v_start, Props(4) - Props(2)*log(-pc0)], 'k--');
% lgd = legend('Simulated CD', 'ICL, $\lambda$ line', 'CSL', 'Swelling $\kappa$ line');
% lgd.FontSize = 14;
% lgd.Interpreter = 'Latex';
% xlabel('$-p''$ (kPa)', 'Interpreter','latex');
% ylabel('Specific volume $\bar{v}$', 'Interpreter','latex');
% set(gca, 'TickLabelInterpreter','Latex');
% set(gca, 'FontSize', 13); % Adjust axis tick number/value size
% set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
% set(get(gca,'YLabel'),'FontSize',15);
% if export
%     export_fig(sprintf('../Cropped_PDF_figures/CDvlogp_OCR%d.png', -pc0/pconf0), '-m7', '-painters');
% end
% 
% 
figure(4);
set(gcf, 'Color', 'w');
plot(-p_data, q_data, 'g', 'LineWidth', 1); grid on; hold on; % Slope = 3
axis equal;
xlabel('$p''$ (kPa)', 'Interpreter','latex');
ylabel('von Mises stress $q = \|\sigma_1 - \sigma_3\|$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);

% 
% 
% function [] = MCC_ys(pc, M)  % Add MCC yield surface plot
%     fplot(@(x) sqrt(x.*(-pc - x))*M, [0, -pc], 'k--'); hold on
% end

