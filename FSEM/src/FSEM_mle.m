function [paramsA, paramsN] = FSEM_mle(y, X, famid, zygosity)
% FSEM_mle does estimation at all tract points
%   using maximum likelihood approach
%
% argin
% y:        N x M matrix
% X:        N x p matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector
%
% argout
% paramsA: (3+p) x M matrix - estimates under the alternative
%          1. The first 3 x M matrix corresponds to 
%          log(sigma_a^2(s)), log(sigma_c^2(s)), and log(sigma_e^2(s))
%          2. The last p x M matrix corresponds to varying coefficients
% paramsN: (2+p) x M matrix - estimated under the null
%          1. The first 2 x M matrix corresponds to 
%          log(sigma_c^2(s)), and log(sigma_e^2(s))
%          2. The last p x M matrix corresponds to varying coefficients
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA


M = size(y, 2);
p = size(X, 2);

options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton');

if nargout == 1
    paramsA = zeros(p+3, M);
    % do analysis at each tract point
    for kk = 1:M
        paramsA(:,kk) = fminunc(@(params) FSEM_lla(params, y(:,kk), X, famid, zygosity), ...
            zeros(p+3,1), options);
    end
end

if nargout == 2
    paramsA = zeros(p+3, M);
    paramsN = zeros(p+2, M);
    % do analysis at each tract point
    for kk = 1:M
        paramsA(:,kk) = fminunc(@(params) FSEM_lla(params, y(:,kk), X, famid, zygosity), ...
            zeros(p+3,1), options);
        paramsN(:,kk) = fminunc(@(params) FSEM_lln(params, y(:,kk), X, famid, zygosity), ...
            zeros(p+2,1), options);
    end
end

end