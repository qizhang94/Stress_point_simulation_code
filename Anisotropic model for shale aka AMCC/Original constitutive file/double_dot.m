function C = double_dot(A,B)
% double contraction between tensors
m = size(A); m(end-1:end) = [];
n = size(B); n(1:2) = [];

if length(m) == 2 && length(n) == 2
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    C(i, j, k, l) = trace(reshape(A(i, j, :, :), [3, 3]) * B(:, :, k, l)');
                end
            end
        end
    end
elseif length(m) == 2 && isempty(n)
    for i = 1:3
        for j = 1:3
            C(i, j) = trace(reshape(A(i, j, :, :), [3, 3]) * B');
        end
    end
elseif length(n) == 2 && isempty(m)
    for i = 1:3
        for j = 1:3
            C(i, j) = trace(A* transpose(reshape(B(:, :, i, j), [3, 3])));
        end
    end
elseif isempty(m) && isempty(n)
    C = trace(A*B');
end                  

end
