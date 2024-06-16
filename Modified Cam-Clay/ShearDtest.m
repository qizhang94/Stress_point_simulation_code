clear; clearvars -global; clc; close all;
restoredefaultpath;

% Investigate the trend of shear (drained) test: shear-induced compaction
Props = [1.2468, 0.2508, 0.05387, 3.5804, 0.25]; % Material properties (v1 = 3.5804 is the value of p = 1 unit, 0.25 is the Poisson's ratio)

% ICL: Props(4) - Props(2)*ln(p/1[unit])
% CSL: Props(4) - (Props(2)-Props(3))*ln(2) - Props(2)*ln(p/1[unit])

% Preconsolidation pressure = 800 [kPa]
pconf0 = 100; % Initial confining pressure (isotropic) [kPa]
pc0 = -800; % [kPa]
stress_old = [-pconf0; -pconf0; -pconf0; 0; 0; 0];
stress_new = stress_old;
v_start = Props(4) - Props(2)*log(-pc0) + Props(3)*log(-pc0/pconf0); % Starting specific volume
hsv_old = [pc0; 0; v_start]; % History variables
hsv_new = hsv_old;
strain_accu = zeros(6,1); % Strain change "from" the isotropic compression state

shear_incr = 2e-3; % Shear strain amount increment (engineering strain)
ntolstep = 500; % Number of total simulation steps

b = [-pconf0; -pconf0; -pconf0; 0; shear_incr; 0]; % RHS of the non-linear equation
q_data = [0]; es_data = [0]; ev_data = [0]; v_data = [v_start]; p_data = [pconf0]; % Resulting data for plot

for i = 1:1:ntolstep
    strain_iter = [0; 0; 0; 0; shear_incr; 0];
    error_tol = 1e-12;
    k = 0; % Outer iteration counter, inner is in the UMAT function
    while 1 > 0
        [stress_new, hsv_new, DDSDDE, DELAS]=MCC_UMAT(Props, stress_old, strain_iter, hsv_old);
        R = stress_new - b; % Residual
        R(5) = 0; % Because the increment should be 0
        if k == 0
            R_ini = R + 1; % Initial residual for this step, the addition of 1 is arbitrary
        end

        if norm(R)/norm(R_ini) < error_tol || norm(R) < error_tol
            break;
        end
        if k > 15  % Not converge
            break;
        end

        DDSDDE(5,:) = 0; DDSDDE(:, 5) = 0; DDSDDE(5,5) = 1;
        strain_iter = strain_iter - DDSDDE\R; % Guarantee the third component is unchanged
        k = k + 1;
    end
    % Save state
    stress_old = stress_new;
    hsv_old = hsv_new;
    strain_accu = strain_accu + strain_iter;

    q_data = [q_data; sqrt(3)*abs(stress_old(5))];
    es_data = [es_data; strain_accu(5)];
    ev_data = [ev_data; sum(strain_accu(1:3))];
    v_data = [v_data; hsv_old(3)];
    p_data = [p_data; -sum(stress_old(1:3))/3]; % not changed

end

figure(1);
set(gcf, 'Color', 'w');
plot(es_data,  q_data, 'b', 'LineWidth', 1); grid on;
xlabel('$\gamma$', 'Interpreter','latex'); 
ylabel('von Mises stress $q = \sqrt{3}\|\sigma_{xz}\|$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


figure(2); % when OCR = 8, shear-induced dilation (epsilon_v > 0)
set(gcf, 'Color', 'w');
plot(es_data, ev_data, 'r', 'LineWidth', 1); grid on;
% set(gca, 'YDir','reverse');
xlabel('$\gamma$', 'Interpreter','latex');
ylabel('$\varepsilon_v$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


figure(3);
set(gcf, 'Color', 'w');
semilogx(p_data, v_data, 'g', 'LineWidth', 1); grid on; hold on;
semilogx([p_data; -pc0], Props(4) - Props(2)*log([p_data; -pc0]), 'b-.', 'LineWidth', 1);
semilogx([p_data; -pc0], Props(4)  - (Props(2)-Props(3))*log(2) - Props(2)*log([p_data; -pc0]), 'r-.', 'LineWidth', 1);
semilogx([pconf0, -pc0],[v_start, Props(4) - Props(2)*log(-pc0)], 'k--');
lgd = legend('Simulated shear', 'ICL, $\lambda$ line', 'CSL', 'Swelling $\kappa$ line', 'Interpreter', 'Latex');
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
xlabel('$p''$ (kPa)', 'Interpreter','latex');
ylabel('Specific volume', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


figure(4);
set(gcf, 'Color', 'w');
plot(p_data, q_data, 'g', 'LineWidth', 1); grid on; hold on; % Slope = 3
plot([0; p_data], Props(1)*[0; p_data], 'r', 'LineWidth', 1);
MCC_ys(pc0, Props(1)); MCC_ys(hsv_new(1), Props(1));
lgd = legend('Simulated shear', 'CSL', 'Yield surface', 'Location','best', 'Interpreter', 'Latex'); axis equal;
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
xlabel('$p''$ (kPa)', 'Interpreter','latex');
ylabel('von Mises stress $q = \sqrt{3}\|\sigma_{xz}\|$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


function [] = MCC_ys(pc, M)  % Add MCC yield surface plot
    fplot(@(x) sqrt(x.*(-pc - x))*M, [0, -pc], 'k--'); hold on
end

