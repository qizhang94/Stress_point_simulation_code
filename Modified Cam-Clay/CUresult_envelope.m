clear; clearvars -global; clc; close all;
restoredefaultpath;
addpath('C:\export_fig-master');

% Effective stress data (NOTE that the confining stress also changes!)
D1 = [1.730250808262685e+02, 5.424600666355130e+02]; % dataset 1: sigma_3 < sigma_1
D2 = [2.549591593869848e+02, 7.988960037038393e+02];
D3 = [2.710978580464739e+02, 8.494618380951232e+02];

% Pore water pressure
PWP1 = -73.025080826268493; % negative pore pressure when OCR is large
PWP2 = 3.450408406130152e+02;
PWP3 = 5.289021419535261e+02;

[c, phi] = ff((D1(1) + D1(2))/2, 0, (D1(2) - D1(1))/2, (D2(1) + D2(2))/2, 0, (D2(2) - D2(1))/2);

figure(1);
set(gcf, 'Color', 'w');
h1 = viscircles([(D1(1) + D1(2))/2, 0], (D1(2) - D1(1))/2, 'Color', 'r'); hold on; grid on;
h2 = viscircles([(D2(1) + D2(2))/2, 0], (D2(2) - D2(1))/2, 'Color', 'y');
h3 = viscircles([(D3(1) + D3(2))/2, 0], (D3(2) - D3(1))/2, 'Color', 'g'); axis equal
xlim([0, inf]);
lgd = legend([h1 h2 h3],{'OCR = 8','OCR = 4/3','OCR = 1'});
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
lgd.Location = 'best';
xlabel('Normal compressive stress $\sigma_n$ (kPa)', 'Interpreter','latex'); 
ylabel('Shear stress $\tau$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
box on;
title(['$c'' =$ ', num2str(c), ' (kPa) and $\phi'' = $', num2str(phi)], 'fontsize', 18, 'Interpreter','latex');
% export_fig ../Cropped_PDF_figures/CU_envelope1.png -m7 -painters;

[c, phi] = ff((D1(1) + D1(2))/2 + PWP1, 0, (D1(2) - D1(1))/2, (D3(1) + D3(2))/2 + PWP3, 0, (D3(2) - D3(1))/2);
figure(2);
set(gcf, 'Color', 'w');
h1 = viscircles([(D1(1) + D1(2))/2 + PWP1, 0], (D1(2) - D1(1))/2, 'Color', 'r'); hold on; grid on;
h2 = viscircles([(D2(1) + D2(2))/2 + PWP2, 0], (D2(2) - D2(1))/2, 'Color', 'y');
h3 = viscircles([(D3(1) + D3(2))/2 + PWP3, 0], (D3(2) - D3(1))/2, 'Color', 'g'); axis equal
xlim([0, inf]);
lgd = legend([h1 h2 h3],{'OCR = 8','OCR = 4/3','OCR = 1'});
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
lgd.Location = 'best';
xlabel('Normal compressive stress $\sigma_n$ (kPa)', 'Interpreter','latex'); 
ylabel('Shear stress $\tau$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
box on;
title(['$c_{\rm u} =$ ', num2str(c), ' (kPa) and $\phi_{\rm u} = $', num2str(phi)], 'fontsize', 18, 'Interpreter','latex');
% export_fig ../Cropped_PDF_figures/CU_envelope2.png -m7 -painters;

% By using this approach, for the effective Mohr's circle, the phi value is
% correct, but c is close to 0. For the total Mohr's circle, the phi value
% is small, but the c is large

function [c, phi] = ff(x1, y1, r1, x2, y2, r2)
    Deltaplus = (x1 - x2)^2 + (y1 - y2)^2 - (r1 + r2)^2;
    Deltaminus = (x1 - x2)^2 + (y1 - y2)^2 - (r1 - r2)^2;
    p1 = r1*(x2^2 + y2^2 - x1*x2 - y1*y2);
    p2 = r2*(x1^2 + y1^2 - x1*x2 - y1*y2);
    q = x1*y2 - x2*y1;
    A = (x2 - x1)*(r1 - r2) + (y1 - y2)*sqrt(Deltaminus);
    B = (y2 - y1)*(r1 - r2) + (x2 - x1)*sqrt(Deltaminus);
    C = -p1 - p2 + q*sqrt(Deltaminus);
    
    c = -C/B;
    phi = atan(-A/B); % in radian
end