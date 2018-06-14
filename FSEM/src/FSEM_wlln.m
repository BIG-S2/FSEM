function [f] = FSEM_wlln(log_sigmaN, R, famid, zygosity, sM, s, h)
% FSEM_wlln_a computes the - weighted log likelihood 
%   under the null (sigma_a2 = 0) at given tract point
% 
% argin
% log_sigmaN: 2 x 1 vector
% R:          N x M matrix <- Residuals: R = y - X*beta
% famid:      N x 1 vector
% zygosity:   N x 1 vector
% sM:         1 x M vector; location points
% s:          scalar; location point
% h:          scalar; bandwidth
%
% argout
% f: - weighted log likelihood under null
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA


% transformation so that sigma's are non-negative
sigma_c2 = exp(log_sigmaN(1));
sigma_e2 = exp(log_sigmaN(2));

a = sigma_c2+sigma_e2; 
b = sigma_c2;

[MZtp1, MZtp2, DZtp1, DZtp2, MDti] = FSEM_index(famid, zygosity);

n1 = sum(MZtp1); n2 = sum(DZtp1); n3 = sum(MDti);

M1 = sum(R(MZtp1,:).^2 + R(MZtp2,:).^2) + sum(R(DZtp1,:).^2 + R(DZtp2,:).^2);
M2 = sum(R(MZtp1,:).*R(MZtp2,:)) + sum(R(DZtp1,:).*R(DZtp2,:));
M0 = sum(R(MDti,:).^2);

% kernel weight vector 1 x M
w = ker((sM-s)/h)/h; sw = sum(w); 
M = length(sM);

f = ( (a*dot(M1,w)-2*b*dot(M2,w))/(a^2-b^2) + (n1+n2)*log(a^2-b^2)*sw + ...
      dot(M0,w)/a + n3*log(a)*sw )/2/M;

end