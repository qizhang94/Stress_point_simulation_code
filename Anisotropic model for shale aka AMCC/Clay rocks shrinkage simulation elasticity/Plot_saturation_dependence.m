%% Opalinus shale
clear; clearvars -global; clc; close all;

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


figure(1); set(gcf, 'Color', 'w');
plot(Eh_expr(:,1), Eh_expr(:,2), 'ro', 'LineWidth',1, 'MarkerSize', 6, 'MarkerFaceColor', 'r'); hold on;
plot(Ev_expr(:,1), Ev_expr(:,2), 'bd', 'LineWidth',1, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
plot(Gvh_expr(:,1), Gvh_expr(:,2), 'k*', 'LineWidth',1, 'MarkerSize', 7.5, 'MarkerFaceColor', 'k');


xlabel('Degree of saturation', 'Interpreter','latex');
ylabel('Stiffness moduli (GPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',16); % Adjust label size
set(get(gca,'YLabel'),'FontSize',16);


mdl_Eh = fitlm(Eh_expr(:,1), Eh_expr(:,2));
mdl_Ev = fitlm(Ev_expr(:,1), Ev_expr(:,2));
mdl_Gvh = fitlm(Gvh_expr(:,1), Gvh_expr(:,2));

plot(Eh_expr(:,1), mdl_Eh.Fitted, 'r--', 'LineWidth',1.2);
plot(Ev_expr(:,1), mdl_Ev.Fitted, 'b--', 'LineWidth',1.2);
plot(Gvh_expr(:,1), mdl_Gvh.Fitted, 'k--', 'LineWidth',1.2);


figure(2); set(gcf, 'Color', 'w');
plot(nuhh_expr(:,1), nuhh_expr(:,2), 'ro', 'LineWidth',1, 'MarkerSize', 6, 'MarkerFaceColor', 'r'); hold on;
plot(nuvh_expr(:,1), nuvh_expr(:,2), 'bd', 'LineWidth',1, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
xlabel('Degree of saturation', 'Interpreter','latex');
ylabel('Poisson''s ratio', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',16); % Adjust label size
set(get(gca,'YLabel'),'FontSize',16);

mdl_nuhh = fitlm(nuhh_expr(:,1), nuhh_expr(:,2));
mdl_nuvh = fitlm(nuvh_expr(:,1), nuvh_expr(:,2));

plot(nuhh_expr(:,1), mdl_nuhh.Fitted, 'r--', 'LineWidth',1.2);
plot(nuvh_expr(:,1), mdl_nuvh.Fitted, 'b--', 'LineWidth',1.2);