
function [LRT_vec, p_vec] = FSEM_lrt(paramsA, paramsN, y, X, famid, zygosity)
% FSEM_lrt does pointwise test given parameter estimates via LRT
%
% argin
% paramsA: (3+p) x M matrix <- output given by FSEM_mle
% paramsN: (2+p) x M matrix <- output given by FSEM_mle
% y:        N x M matrix
% X:        N x p matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector
%
% argout
% LRT_vec: 1 x M vector; test statistics
% p_vec:   1 x M vector; p values
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA


K = size(paramsA, 2);

LRT_vec = zeros(1, K);

for kk = 1:K
    fa = FSEM_lla(paramsA(:,kk), y(:,kk), X, famid, zygosity);
    fn = FSEM_lln(paramsN(:,kk), y(:,kk), X, famid, zygosity);
    LRT_vec(kk) = 2*max(0, fn-fa);
end

p_vec = 1-normcdf(sqrt(LRT_vec));

end
