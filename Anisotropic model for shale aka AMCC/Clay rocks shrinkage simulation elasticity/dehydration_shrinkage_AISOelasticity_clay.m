%% Opalinus shale
clear; clearvars -global; clc; close all;

% Saturation = 0.44; nu_hh = 0.1744, nu_vh = 0.189638

Ev_sat = 190; % MPa

% Constant stiffness anisotropy uses constant elastic parameters reported for Opalinus shale at 44\% saturation
% Biot tensor: [0.56, 0.56, 0.77, 0, 0, 0]' for constant stiffness, corrsponding Ks is 41 [GPa]
% Biot tensor: [0.49, 0.49, 0.8, 0, 0, 0]' for variable stiffnes at 10\% saturation, corrsponding Ks is 38 [GPa]

% Ev = linf(Ev_sat,0.44,4595); % MPa
% Eh = linf(Eh_sat,0.44,23069); % MPa
% nuhh = 0.1744;
% nuvh = 0.189638;
% Gvh = linf(Gvh_sat,0.44,2780); % MPa

oneVoigt = [1,1,1,0,0,0]';
Ks = inf; % Solid grain bulk modulus
figure(1); set(gcf, 'Color', 'w');
for suction_stress = linspace(5,105,20) % MPa
    psi_w = SWCC(2, 0, 1, 1/16.67, suction_stress); % Obtain water saturation
    Ev = linf(Ev_sat,psi_w,-60); % MPa
    Eh = 300; % MPa
    nuhh = 0.125;
    nuvh = 0.125;
    Gvh = 150; % MPa

% It seems that picewise relation does not perform better than constant stiffness

%     Ev = interp1(Ev_expr(:,1), Ev_expr(:,2), psi_w, "linear", "extrap")*1e3;
%     Eh = interp1(Eh_expr(:,1), Eh_expr(:,2), psi_w, "linear", "extrap")*1e3;
%     Gvh = interp1(Gvh_expr(:,1), Gvh_expr(:,2), psi_w, "linear", "extrap")*1e3;
%     nuhh = interp1(nuhh_expr(:,1), nuhh_expr(:,2), psi_w, "linear", "extrap");
%     nuvh = interp1(nuvh_expr(:,1), nuvh_expr(:,2), psi_w, "linear", "extrap");
    
    
    % Axis of cross anisotropy coincides with the z-axis
    S11 = [1/Eh, -nuhh/Eh, -nuvh/Ev; ...
        -nuhh/Eh, 1/Eh, -nuvh/Ev; ...
        -nuvh/Ev, -nuvh/Ev, 1/Ev];
    S = [S11,zeros(3,3);zeros(3,3),diag([2*(1+nuhh)/Eh, 1/Gvh, 1/Gvh])]; % Compliance matrix 6*6 (VTI)
    Ce = S\eye(6); % Stiffness matrix 6*6 (VTI), MPa
        
    alphaBiot = oneVoigt; % Voigt notation: 6*1 vector
    strain_vec = S*(psi_w*-suction_stress*alphaBiot); % compressive strain has the negative sign
    
    % However, when plotting, we adopt the sign convention of soil mechanics
    %semilogx(suction_stress, -strain_vec(1), 'bo'); grid on; hold on; % BP = bed-parallel
    %semilogx(suction_stress, -strain_vec(3), 'rs'); % BN = bed-normal
    semilogx(suction_stress, strain_vec(3)/strain_vec(1), 'bd'); grid on; hold on; % BN = bed-normal
end


% Soil Water Characteristic Curve (Water Retention Law)

function psi_w = SWCC(n, psi_res, psi_max, alpha, suction)
    m = 1-1/n;
    psi_w = psi_res + (psi_max - psi_res)./(1+(alpha*suction).^n).^m;
end

function Xpsiw = linf(Xsat, psi_w, aX)
    Xpsiw = Xsat + (1-psi_w)*aX;
    % aX has the same unit as Xsat and Xpsiw
end