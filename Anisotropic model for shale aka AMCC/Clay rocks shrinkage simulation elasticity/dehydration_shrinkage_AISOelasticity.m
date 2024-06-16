%% Opalinus shale
clear; clearvars -global; clc; close all;

% Saturation = 0.44; nu_hh = 0.1744, nu_vh = 0.189638

Ev_sat = 5222; % MPa
Eh_sat = 13615; % MPa
Gvh_sat = 3559; % MPa
nuvh_sat = 0.395;
nuhh_sat = -0.26;

Ev_expr = [0.05625250866630175, 9.233716475095783
0.23919722678343372, 8.314176245210724
0.44494070425104904, 8.007662835249034
0.9442218573253055, 5.555555555555543]; % GPa

Eh_expr = [0.05359605911330055, 37.12643678160919
0.23718299580368557, 29.463601532567047
0.435613939062215, 25.938697318007662
0.94519978106185, 15.287356321839086];

Gvh_expr = [0.05846378398102534, 6.015325670498084
0.24140120415982483, 5.172413793103444
0.44334245575624887, 4.7892720306513406
0.9406166757890895, 3.4099616858237436];

nuhh_expr = [0.05534351145038173, 0.16644295302013423
0.24236641221374036, 0.19865771812080532
0.4408396946564886, 0.17583892617449667
0.9561068702290078, -0.29798657718120797];

nuvh_expr = [0.05725190839694652, 0.12483221476510066
0.24045801526717553, 0.15436241610738255
0.4408396946564886, 0.18926174496644294
0.952290076335878, 0.4];

% Constant stiffness anisotropy uses constant elastic parameters reported for Opalinus shale at 44\% saturation
% Biot tensor: [0.56, 0.56, 0.77, 0, 0, 0]' for constant stiffness, corrsponding Ks is 41 [GPa]
% Biot tensor: [0.49, 0.49, 0.8, 0, 0, 0]' for variable stiffnes at 10\% saturation, corrsponding Ks is 38 [GPa]

% Ev = linf(Ev_sat,0.44,4595); % MPa
% Eh = linf(Eh_sat,0.44,23069); % MPa
% nuhh = 0.1744;
% nuvh = 0.189638;
% Gvh = linf(Gvh_sat,0.44,2780); % MPa

oneVoigt = [1,1,1,0,0,0]';
Ks = 41E3; % Solid grain bulk modulus
figure(1); set(gcf, 'Color', 'w');
sigma_total_conf = [0; 0; 0; 0; 0; 0]; % 2 [MPa] confining stress
for suction_stress = linspace(0,400,51) % MPa
    psi_w = SWCC(2.083, 0, 1, 1/47.4, suction_stress); % Obtain water saturation
    Ev = linf(Ev_sat,psi_w,4595); % MPa
    Eh = linf(Eh_sat,psi_w,23069); % MPa
    nuhh = linf(nuhh_sat,psi_w,0.557);
    nuvh = linf(nuvh_sat,psi_w,-0.31);
    Gvh = linf(Gvh_sat,psi_w,2780); % MPa

% Ev = linf(Ev_sat,0.44,4595); % MPa
% Eh = linf(Eh_sat,0.44,23069); % MPa
% nuhh = 0.1744;
% nuvh = 0.189638;
% Gvh = linf(Gvh_sat,0.44,2780); % MPa

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
        
    alphaBiot = oneVoigt - Ce*oneVoigt/3/Ks; % Voigt notation: 6*1 vector
    strain_vec = S*(sigma_total_conf + psi_w*-suction_stress*alphaBiot); % compressive strain has the negative sign
    
    % However, when plotting, we adopt the sign convention of soil mechanics
    plot(suction_stress, strain_vec(1), 'bo', 'LineWidth', 1,  'MarkerFaceColor','w'); grid off; hold on; % BP = bed-parallel
    plot(suction_stress, strain_vec(3), 'rs', 'LineWidth', 1, 'MarkerSize', 7.5, 'MarkerFaceColor','y'); % BN = bed-normal
end


xlabel('Suction (MPa)', 'Interpreter','latex');
ylabel('Strain', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',16); % Adjust label size
set(get(gca,'YLabel'),'FontSize',16);
lgd = legend('BP direction','BN direction', 'location', 'best');
legend boxoff
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
title('Variable stiffness','FontSize',16, 'Interpreter','latex');

% Soil Water Characteristic Curve (Water Retention Law)

function psi_w = SWCC(n, psi_res, psi_max, alpha, suction)
    m = 1-1/n;
    psi_w = psi_res + (psi_max - psi_res)./(1+(alpha*suction).^n).^m;
end

function Xpsiw = linf(Xsat, psi_w, aX)
    Xpsiw = Xsat + (1-psi_w)*aX;
    % aX has the same unit as Xsat and Xpsiw
end