clear; clearvars -global; clc; close all;
restoredefaultpath;
addpath('C:\export_fig-master');

% Effective stress data
D1 = [100, 3.134493788186313e+02]; % dataset 1: sigma_3 < sigma_1
D2 = [600, 1.878009487245854e+03];
D3 = [800, 2.502705217251116e+03];

% Shear data
DS1 = [27.957380933430414, 1.720426190665696e+02];
DS2 = [1.683048062769794e+02, 1.031695193723021e+03];

% Pore water pressure
PWP1 = 0;
PWP2 = 0;
PWP3 = 0;

[c, phi] = ff((D1(1) + D1(2))/2 + PWP1, 0, (D1(2) - D1(1))/2, (D2(1) + D2(2))/2 + PWP2, 0, (D2(2) - D2(1))/2);
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
% export_fig ../Cropped_PDF_figures/CD_envelope1.png -m7 -painters;

[c, phi] = ff((DS1(1) + DS1(2))/2, 0, (DS1(2) - DS1(1))/2, (DS2(1) + DS2(2))/2, 0, (DS2(2) - DS2(1))/2);
figure(2);
set(gcf, 'Color', 'w');
h1 = viscircles([(DS1(1) + DS1(2))/2, 0], (DS1(2) - DS1(1))/2, 'Color', 'b'); hold on; grid on;
h2 = viscircles([(DS2(1) + DS2(2))/2, 0], (DS2(2) - DS2(1))/2, 'Color', 'm'); axis equal
xlabel('Normal compressive stress $\sigma_n$ (kPa)', 'Interpreter','latex'); 
ylabel('Shear stress $\tau$ (kPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);
lgd = legend([h1 h2],{'OCR = 8','OCR = 4/3'});
lgd.FontSize = 14;
lgd.Interpreter = 'Latex';
lgd.Location = 'best';
box on;
title(['$c'' =$ ', num2str(c), ' (kPa) and $\phi'' = $', num2str(phi)], 'fontsize', 18, 'Interpreter','latex');


% By using this approach, we CAN get a correct phi value (triaxial drained), but c = 0.
% 通过这种方法，得到的包络线与CSE578-Lecture-2-Example-2.2中的MC屈服面并不一致：斜率相同，但没有截距

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