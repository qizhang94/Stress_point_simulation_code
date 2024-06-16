function detA = detect_bifur(Cep)
% Cep is a 9*9 matrix, equivalent to a fourth order continuum elastoplastic
% tensor. We assume the normal vector of the shear band could be
% represented as n = [cos(psi), 0, sin(psi)]^T

% A^ep_jk = ni*C^ep_ijkl*nl

detA = [];
for psi = linspace(0, pi, 10001)
    Aep = zeros(3,3);
    n = [cos(psi); 0; sin(psi)];
    for j = 1:3
        for k = 1:3
            for i = 1:3
                for l = 1:3
                    Cijkl = (Cep((j-1)*3+i,(l-1)*3+k) + Cep((i-1)*3+j,(l-1)*3+k)...
                        + Cep((j-1)*3+i,(k-1)*3+l) + Cep((i-1)*3+j,(k-1)*3+l))/4; % ensure minor-symmetry
                    Aep(j,k) = Aep(j,k) + n(i)*Cijkl*n(l);
                end
            end
        end
    end
    detA = [detA;nthroot(det(Aep),3)];
end


end