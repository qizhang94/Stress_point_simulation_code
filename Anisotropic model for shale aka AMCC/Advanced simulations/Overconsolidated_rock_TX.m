clear; clearvars -global; clc; close all;
export = false;
% if export
%     addpath('C:\export_fig-master');
% end

% Investigate the trend of triaxial drained test
Props = [1.07, 0.0026, 0, pi/6]; % M, lambda, kappa, theta


% Preconsolidation pressure = 40 [MPa]
pconf0 = 1; % Initial confining pressure (isotropic) [MPa]
pc0 = -40; % [MPa]
stress_old = [-pconf0; -pconf0; -pconf0; 0; 0; 0];
stress_new = stress_old;
hsv_old = [pc0; 0]; % History variables
hsv_new = hsv_old;
strain_accu = zeros(6,1); % Strain change "from" the isotropic compression state

comp_incr = -1e-4; % Compression strain amount increment
ntolstep = 500; % Number of total simulation steps

b = [-pconf0; -pconf0; comp_incr; 0; 0; 0]; % RHS of the non-linear equation
q_data = 0; ea_data = 0; ev_data = 0; p_data = pconf0; % Resulting data for plot

flag = false;
figure;
set(gcf, 'Color', 'w');
for i = 1:1:ntolstep
    strain_iter = [0; 0; comp_incr; 0; 0; 0];
    error_tol = 1e-12;
    k = 0; % Outer iteration counter, inner is in the UMAT function
    while 1 > 0
        [stress_new, hsv_new, DDSDDE, DELAS, Cep]=AMCC_UMAT_bifurcation(Props, stress_old, strain_iter, hsv_old); % Cep: 9*9
        R = stress_new - b; % Residual
        R(3) = 0; % Because we know the epsilon_zz, then the increment should be 0
        if k == 0
            R_ini = R + 1; % Initial residual for this step, the addition of 1 is arbitrary
        end

        if norm(R)/norm(R_ini) < error_tol || norm(R) < error_tol
            break;
        end
        if k > 15  % Not converge
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

    % Detect bifurcation
    detA = detect_bifur(Cep);
    psideg = linspace(0, 180, length(detA)); % unit of degree
    [~,I] = min(detA); % return the position to get psi
    disp(['e_a = ' num2str(-strain_accu(3)), ': det = ',  num2str(min(detA)),...
        ', psi = ', num2str(min(psideg(I), 180-psideg(I)))]);
    if min(detA) < 0 && flag == false
        flag = true;
        semilogy(psideg, 1e3+(detA)); xlabel('Psi, degree'); ylabel('det(A^ep)'); % an offset of 10^3
    end

end

figure;
set(gcf, 'Color', 'w');
plot(ea_data,  q_data, 'b', 'LineWidth', 1); grid on;
xlabel('Compressive axial strain $\varepsilon_a$', 'Interpreter','latex'); 
ylabel('von Mises stress $q = \|\sigma_1 - \sigma_3\|$ (MPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
% if export
%     export_fig(sprintf('../Cropped_PDF_figures/CDeaq_OCR%d.png', -pc0/pconf0), '-m7', '-painters');
% end

figure;
set(gcf, 'Color', 'w');
plot(ea_data, -ev_data, 'r', 'LineWidth', 1); grid on;
% set(gca, 'YDir','reverse');
xlabel('Compressive axial strain $\varepsilon_a$', 'Interpreter','latex');
ylabel('Volumetric strain $\varepsilon_v$', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
% if export
%     export_fig(sprintf('../Cropped_PDF_figures/CDeaev_OCR%d.png', -pc0/pconf0), '-m7', '-painters');
% end

