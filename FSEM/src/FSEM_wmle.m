function [betaA, paramsA, paramsN] = ...
    FSEM_wmle(y, X, famid, zygosity, sM, h)
% FSEM_wmle does estimation at all tract points
%   using weighted maximum likelihood approach
%
% argin
% y:        N x M matrix
% X:        N x p matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector
% sM:       1 x M vector; tract points
% h:        scalar; bandwidth
%
% argout
% betaA:   p x M matrix; varying coefficients matrix
% paramsA: 3 x M matrix - estimates under the alternative
% paramsN: 2 x M matrix - estimated under the null
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA


loc = sM;

if nargin < 6
    h = 5*range(sM)/length(sM)*ones(1,length(loc));
end

if length(h) == 1
    h = h*ones(1, length(loc));
end

K = length(loc);

options = optimoptions('fminunc','GradObj','off', 'Display','off','Algorithm','quasi-newton');

paramsA = FSEM_mle(y, X, famid, zygosity);
betaA = paramsA(4:end,:);
R = y - X*betaA;

if nargout == 2
    paramsA = zeros(3, K);
    for kk = 1:K
        paramsA(:, kk) = fminunc(@(params) ...
            FSEM_wlla(params, R, famid, zygosity, sM, loc(kk), h(kk)), ...
            zeros(3,1), options);
    end
end

if nargout == 3
    paramsA = zeros(3, K);
    paramsN = zeros(2, K);
    for kk = 1:K
        paramsA(:,kk)  = fminunc(@(params) ...
            FSEM_wlla(params, R, famid, zygosity, sM, loc(kk), h(kk)), ...
            zeros(3,1), options);
        paramsN(:, kk) = fminunc(@(params) ...
            FSEM_wlln(params, R, famid, zygosity, sM, loc(kk), h(kk)), ...
            zeros(2,1), options);
    end
end

end