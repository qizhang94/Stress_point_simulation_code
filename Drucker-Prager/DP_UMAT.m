function [STRESS, hsv, DDSDDE]=DP_UMAT(PROPS, STRESS0, DSTRAIN0, hsv0)
% Return mapping and CTO of the DP plasticity algorithm (2D plane strain)
% STRESS0: Old
% STRESS: New
% DDSDDE: Algorithmic consistent tangent operator
% Assume Ce (6*6) is uniform and constant


E = PROPS(1); % Young's modulus
nu = PROPS(2); % Poisson's ratio
phi = PROPS(3);  %  phi  -  angle of friction (rad)
psi = PROPS(4);  %  psi  -  angle of dilation (rad)
coh = PROPS(5);  %  c  -  cohesion

mu = E/2/(1+nu);
K = E/3/(1-2*nu);

Ce = K*[1;1;1;0;0;0]*[1,1,1,0,0,0] + 2*mu*(diag([1,1,1,0.5,0.5,0.5]) - ...
    1/3*[1;1;1;0;0;0]*[1,1,1,0,0,0]);

% A = 3*sqrt(2)*coh*cos(phi)/sqrt(9+3*(sin(phi))^2);
% B = 3*sqrt(2)*sin(phi)/sqrt(9+3*(sin(phi))^2);
% b = 3*sqrt(2)*sin(psi)/sqrt(9+3*(sin(psi))^2);

A = 2*sqrt(6)*coh*cos(phi)/(3 - sin(phi)); % Pass through the compression corner
B = 2*sqrt(6)*sin(phi)/(3 - sin(phi));
b = 2*sqrt(6)*sin(psi)/(3 - sin(psi));

stress_trial = STRESS0 + Ce*DSTRAIN0;   % Trial stress
stress_trial_tensor = [stress_trial(1), stress_trial(4), stress_trial(5); stress_trial(4), stress_trial(2), stress_trial(6); stress_trial(5), stress_trial(6), stress_trial(3)];
p_trial = 1/3*trace(stress_trial_tensor);
stress_dev_trial = stress_trial_tensor - p_trial*eye(3);
q_trial = sqrt(1.5)*norm(stress_dev_trial, 'fro');

if (sqrt(2/3)*q_trial + B*p_trial - A <= 0)
    % Not yield
    DDSDDE = Ce;
    STRESS = stress_trial;
    hsv = hsv0;
else
    % Yield, more complicated
    one_voigt = [1;1;1;0;0;0];
    I4 = diag([1,1,1,0.5,0.5,0.5]);
    nhat_tensor = stress_dev_trial/norm(stress_dev_trial, 'fro');
    nhat = [nhat_tensor(1,1); nhat_tensor(2,2); nhat_tensor(3,3); nhat_tensor(1,2); nhat_tensor(1,3); nhat_tensor(2,3)];
    delta_plastic = (sqrt(2/3)*q_trial + B*p_trial - A)/(2*mu + B*K*b);
    STRESS = stress_trial - delta_plastic*(K*b*one_voigt + 2*mu*nhat);
    hsv = hsv0 + sqrt(2/3)*delta_plastic; % deviatoric strain! Note that norm(nhat_tensor, 'fro') = 1
    DDSDDE = Ce - (K*b*one_voigt + 2*mu*nhat)*...
        transpose(K*B*one_voigt + 2*mu*nhat)/(2*mu + B*K*b)...
        - 4*mu^2*delta_plastic/norm(stress_dev_trial, 'fro')*(I4 - 1/3*(one_voigt*one_voigt') - nhat*nhat');

    if q_trial - sqrt(6)*mu*delta_plastic < 0
        % Return mapping to the apex must be applied!
        hsv = hsv0;
        STRESS = A/B*one_voigt;
        DDSDDE = zeros(6,6);
    end
end

% See if DDSDDE is always invertible (it seems not!)
% coeff = 1/(sum(DDSDDE(1:3,1:3), 'all')/9);

% From 3D to 2D

% DDSDDE = [DDSDDE(1:2, 1:2), DDSDDE(1:2, 4); DDSDDE(4, 1:2), DDSDDE(4,4)];

end