function B = stiffness_to_mat6by6_AppendixB(A)
% Transform a symmetric (minor and major) fourth order tensor A (3*3*3*3) to a matrix B (6*6)
% Only apply to stiffness!
% Note Voigt order: 11 22 33 12 13 23 is used here

e1 = [1, 0, 0; 0, 0, 0; 0, 0, 0];
e2 = [0, 0, 0; 0, 1, 0; 0, 0, 0];
e3 = [0, 0, 0; 0, 0, 0; 0, 0, 1];
e4 = [0, 0.5, 0; 0.5, 0, 0; 0, 0, 0];
e5 = [0, 0, 0.5; 0, 0, 0; 0.5, 0, 0];
e6 = [0, 0, 0; 0, 0, 0.5; 0, 0.5, 0];

D_aug = [reshape(double_dot(A, e1), [9, 1]), reshape(double_dot(A, e2), [9, 1]), ...
    reshape(double_dot(A, e3), [9, 1]), reshape(double_dot(A, e4), [9, 1]),...
    reshape(double_dot(A, e5), [9, 1]), reshape(double_dot(A, e6), [9, 1])]; % 9*6 matrix

B = [D_aug(1,:); D_aug(5,:); D_aug(9,:); D_aug(4,:); D_aug(7,:); D_aug(8,:)];
end

