function [f] = FSEM_lln(params, y, X, famid, zygosity)
% FSEM_lln computes -loglikelihood 
%   under null (sigma_a2 = 0) at given tract point
%
% argin
% paramsA: (2+p) x 1 vector
% y:        N x 1 vector 
% X:        N x p matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector 
%
% argout
% f: -loglikelihood under null

% The first two elements correspond to the variance components,
% transformation so that sigma's are non-negative
sigma_c2 = exp(params(1));
sigma_e2 = exp(params(2));
% The next p elements correspond to covariates' coefficients
beta = params(3:end);

R = y-X*beta;

a = sigma_c2+sigma_e2; 
b = sigma_c2;

[MZtp1, MZtp2, DZtp1, DZtp2, MDti] = FSEM_index(famid, zygosity);

n1 = sum(MZtp1); n2 = sum(DZtp1); n3 = sum(MDti);

SR = sum(R(MZtp1).^2 + R(MZtp2).^2) + sum(R(DZtp1).^2 + R(DZtp2).^2);
XR = 2*sum(R(MZtp1).*R(MZtp2)) + 2*sum(R(DZtp1).*R(DZtp2));
SMD = sum(R(MDti).^2);

f = (a*SR-b*XR)/(a^2-b^2) + (n1+n2)*log(a^2-b^2) + SMD/a + n3*log(a);

f = f/2.0;

end