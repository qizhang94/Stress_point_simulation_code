clear; clearvars -global; clc; close all;


E_hard = 21700; % MPa
nu_hard = 0.23;

E_soft = 11900; % MPa
nu_soft = 0.21;

T_hard = layer(E_hard, nu_hard);
T_soft = layer(E_soft, nu_soft);

A = (T_hard.t1+T_soft.t1)/2 + (2/(T_hard.t3+T_soft.t3))*((T_hard.t4+T_soft.t4)^2/4);
B = (T_hard.t2+T_soft.t2)/2 + (2/(T_hard.t3+T_soft.t3))*((T_hard.t4+T_soft.t4)^2/4);
CC = 2/(T_hard.t3+T_soft.t3);
F = (2/(T_hard.t3+T_soft.t3))*((T_hard.t4+T_soft.t4)/2);
D = 2/(T_hard.t5+T_soft.t5);
M = (T_hard.t6+T_soft.t6)/2;

C = [A, B, F, 0, 0, 0;
    B, A, F, 0, 0, 0;
    F, F, CC, 0, 0, 0;
    0, 0, 0, M, 0, 0;
    0, 0, 0, 0, D, 0;
    0, 0, 0, 0, 0, D]; % Transversely isotropic


lambda = C(1,2); mu_T = C(4,4); mu_L = C(5,5);
alpha = C(1,3)-C(1,2); beta = C(3,3) - lambda - 4*mu_L + 2*mu_T - 2*alpha; % Same unit as Ev/Eh/Gvh for these 5 constants

clearvars -except lambda mu_T mu_L alpha beta


for theta = 0:5:90
    Props = [1.65, 0.0026, 0, theta*pi/180]; % M, lambda, kappa, theta
    
    
    % Preconsolidation pressure = 40 [MPa]
    pconf0 = 14; % Initial confining pressure (isotropic) [MPa]
    pc0 = -14; % [MPa]
    stress_old = [-pconf0; -pconf0; -pconf0; 0; 0; 0];
    stress_new = stress_old;
    hsv_old = [pc0; 0]; % History variables
    hsv_new = hsv_old;
    strain_accu = zeros(6,1); % Strain change "from" the isotropic compression state
    
    comp_incr = -1e-4; % Compression strain amount increment
    ntolstep = 500; % Number of total simulation steps
    
    b = [-pconf0; -pconf0; comp_incr; 0; 0; 0]; % RHS of the non-linear equation
    q_data = 0; ea_data = 0; ev_data = 0; p_data = pconf0; % Resulting data for plot
    
    for i = 1:1:ntolstep
        strain_iter = [0; 0; comp_incr; 0; 0; 0];
        error_tol = 1e-12;
        k = 0; % Outer iteration counter, inner is in the UMAT function
        while 1 > 0
            [stress_new, hsv_new, DDSDDE, DELAS, Cep]=AMCC_UMAT_Synthetic(Props, stress_old, strain_iter, hsv_old); % Cep: 9*9
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
    
    end
    disp([theta, -stress_old(3)]);
end

%% Show results

Tien_6MPa_data = [0, 73.03030303030303
15.06493506493507, 68.98989898989899
30, 61.111111111111114
44.93506493506494, 55.15151515151515
59.74025974025975, 48.08080808080808
74.80519480519482, 66.16161616161617
90, 83.23232323232324
];

Tien_14MPa_data = [0, 95.05050505050505
15.06493506493507, 88.08080808080808
29.870129870129873, 81.11111111111111
44.93506493506494, 70.20202020202021
60.00000000000001, 58.88888888888889
74.80519480519482, 90
90, 105.25252525252526];

Simu_6MPa_data = [         0   75.2470
    5.0000   73.7049
   10.0000   69.6619
   15.0000   64.4131
   20.0000   59.1656
   25.0000   54.6455
   30.0000   51.1497
   35.0000   48.7357
   40.0000   47.3731
   45.0000   47.0274
   50.0000   47.6983
   55.0000   49.4324
   60.0000   52.3171
   65.0000   56.4456
   70.0000   61.8243
   75.0000   68.1810
   80.0000   74.6849
   85.0000   79.8034
   90.0000   81.7816];

Simu_14MPa_data = [         0   93.0149
    5.0000   91.9863
   10.0000   89.1986
   15.0000   85.3716
   20.0000   81.2933
   25.0000   77.5690
   30.0000   74.5638
   35.0000   72.4596
   40.0000   71.3367
   45.0000   71.2363
   50.0000   72.1973
   55.0000   74.2698
   60.0000   77.5046
   65.0000   81.9106
   70.0000   87.3630
   75.0000   93.4550
   80.0000   99.3404
   85.0000  103.7486
   90.0000  105.4032];

figure;
set(gcf, 'Color', 'w');
plot(Tien_6MPa_data(:,1), Tien_6MPa_data(:,2), 'ro', 'MarkerSize', 7, 'MarkerFaceColor','r'); hold on;
plot(Tien_14MPa_data(:,1), Tien_14MPa_data(:,2), 'bo', 'MarkerSize', 7, 'MarkerFaceColor','b');

plot(Simu_6MPa_data(:,1), Simu_6MPa_data(:,2), 'r--', 'LineWidth', 1);
plot(Simu_14MPa_data(:,1), Simu_14MPa_data(:,2), 'b--', 'LineWidth', 1);


xlabel('Bedding plane orientation $\theta$ (deg)', 'Interpreter','latex'); 
ylabel('Maximum vertical stress $\|\sigma_z\|$ (MPa)', 'Interpreter','latex');
set(gca, 'TickLabelInterpreter','Latex');
set(gca, 'FontSize', 13); % Adjust axis tick number/value size
set(get(gca,'XLabel'),'FontSize',15); % Adjust label size
set(get(gca,'YLabel'),'FontSize',15);


function T = layer(E, nu)

lambda = E*nu/(1+nu)/(1-2*nu);
G = E/2/(1+nu);

T.t1 = 4*G*(lambda+G)/(lambda+2*G);
T.t2 = 2*G*lambda/(lambda+2*G);
T.t3 = 1/(lambda+2*G);
T.t4 = lambda/(lambda+2*G);
T.t5 = 1/G;
T.t6 = G;

end