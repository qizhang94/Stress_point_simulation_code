%% 相当于博士论文Appendix D
clear; clearvars -global; clc; close all;
Props = [1.07, 0.0026, 0, pi/3];


pconf0 = 10; % MPa
stress_old = [-pconf0; -pconf0; -pconf0; 0; 0; 0];
stress_new = stress_old;
stress_ini = stress_old;

pc0 = -40; % [MPa]
hsv_old = [pc0; 0]; % History variables
hsv_new = hsv_old;


comp_incr = -1e-3; % Compression strain amount increment
strain_iter = [0; comp_incr; 0; 0; 0; 0];

for lstep = 1:5
    [stress_new, hsv_new, DDSDDE, DELAS]=AMCC_UMAT(Props, stress_old, strain_iter, hsv_old); % Cep: 9*9
    stress_old = stress_new;
    hsv_old = hsv_new;
end

strain_final = DELAS\stress_ini + lstep*strain_iter; % Note the convention for shear strain ==> engineering shear strain