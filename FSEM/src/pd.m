function [P] = pd(X, eps)
% X has to be symmetric

if nargin < 2
    eps = 0.0;
end

[V,D] = eig(X);
D = real(D);
V = real(V);
P = V*(D+abs(D))*V'/2;
P = (P+P')/2. + eps*eye(size(X,1));

end