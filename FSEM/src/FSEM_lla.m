function [f] = FSEM_lla(paramsA, y, X, famid, zygosity)
% FSEM_lla computes the -loglikelihood 
%   under the alternative (sigma_a2 ~= 0 && sigma_c2 ~= 0) at given tract point
%
% argin
% paramsA: (3+p) x 1 vector
% y:        N x 1 vector 
% X:        N x p matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector 
%
% argout 
% f: -loglikelihood under alternative

% The first three elements correspond to the variance components,
% transformation so that sigma's are non-negative
sigma_a2 = exp(paramsA(1)); 
sigma_c2 = exp(paramsA(2)); 
sigma_e2 = exp(paramsA(3));
% The next p elements correspond to covariates' coefficients
beta = paramsA(4:end);

R = y - X*beta;

a = sigma_a2+sigma_c2+sigma_e2; 
b = sigma_a2+sigma_c2; 
c = 0.5*sigma_a2+sigma_c2;

[MZtp1, MZtp2, DZtp1, DZtp2, MDti] = FSEM_index(famid, zygosity);

n1 = sum(MZtp1); n2 = sum(DZtp1); n3 = sum(MDti);

SM = sum(R(MZtp1).^2 + R(MZtp2).^2); 
XM = 2*sum(R(MZtp1).*R(MZtp2)); 

SD = sum(R(DZtp1).^2 + R(DZtp2).^2);
XD = 2*sum(R(DZtp1).*R(DZtp2));

SMD = sum(R(MDti).^2);

f = (a*SM-b*XM)/(a^2-b^2) + n1*log(a^2-b^2) + ...
    (a*SD-c*XD)/(a^2-c^2) + n2*log(a^2-c^2) + SMD/a + n3*log(a);

f = f/2.0;

end