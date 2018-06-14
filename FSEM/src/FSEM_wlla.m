function [f] = FSEM_wlla(log_sigmaA, R, famid, zygosity, sM, s, h, der)
% FSEM_wlla computes the - weighted log likelihood
%   under the alternative (sigma_a2 ~= 0 && sigma_c2 ~= 0) at given tract point
% 
% argin
% log_sigmaA: 3 x 1 vector
% R:          N x M matrix <- Residuals: R = y - X*beta
% famid:      N x 1 vector
% zygosity:   N x 1 vector
% sM:         1 x M vector; location points
% s:          scalar; location point
% h:          scalar; bandwidth
%
% argout
% f: - weighted log likelihood under alternative
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA

if nargin < 8
    der = 0;
end

% transformation so that sigma's are non-negative
sigma_a2 = exp(log_sigmaA(1));
sigma_c2 = exp(log_sigmaA(2));
sigma_e2 = exp(log_sigmaA(3));

a = sigma_a2 + sigma_c2 + sigma_e2;
b = sigma_a2 + sigma_c2;
c = 0.5*sigma_a2 + sigma_c2;

[MZtp1, MZtp2, DZtp1, DZtp2, MDti] = FSEM_index(famid, zygosity);

n1 = sum(MZtp1); n2 = sum(DZtp1); n3 = sum(MDti);

SM = R(MZtp1,:).^2 + R(MZtp2,:).^2; XM = 2*R(MZtp1,:).*R(MZtp2,:); % n1 x M
SD = R(DZtp1,:).^2 + R(DZtp2,:).^2; XD = 2*R(DZtp1,:).*R(DZtp2,:); % n2 x M
SMD = R(MDti,:).^2; % n3 x M

% kernel weight vector 1 x M
w = ker((sM-s)/h)/h; sw = sum(w); M = length(w);

if der == 0
    
    SM = bsxfun(@times, SM, w); XM = bsxfun(@times, XM, w); % n1 x M
    SD = bsxfun(@times, SD, w); XD = bsxfun(@times, XD, w); % n2 x M
    SMD = bsxfun(@times, SMD, w); % n3 x M
    
    f = ( (a*sum(SM(:))-b*sum(XM(:)))/(a^2-b^2) + n1*log(a^2-b^2)*sw + ...
        (a*sum(SD(:))-c*sum(XD(:)))/(a^2-c^2) + n2*log(a^2-c^2)*sw + ...
        sum(SMD(:))/a + n3*log(a)*sw )/2/M; % - loglikelihood
end

if der == 1
    
    SM = bsxfun(@times, SM, w); XM = bsxfun(@times, XM, w); % n1 x M
    SD = bsxfun(@times, SD, w); XD = bsxfun(@times, XD, w); % n2 x M
    SMD = bsxfun(@times, SMD, w); % n3 x M
    
    C1 = (a^2+b^2)/(a^2-b^2)^2; C2 = 2*a*b/(a^2-b^2)^2;
    C3 = (a^2+c^2)/(a^2-c^2)^2; C4 = 2*a*c/(a^2-c^2)^2;
    
    Fa = [-C1*sum(SM,2) + C2*sum(XM,2) + 2*a/(a^2-b^2)*sw;
        -C3*sum(SD,2) + C4*sum(XD,2) + 2*a/(a^2-c^2)*sw;
        -1/a^2*sum(SMD,2) + 1/a*sw]; % (n1+n2+n3) x 1
    Fb = [C2*sum(SM,2) - C1*sum(XM,2) - 2*b/(a^2-b^2)*sw;
        zeros(n2+n3,1)]; % (n1+n2+n3) x 1
    Fc = [zeros(n1,1);
        C4*sum(SD,2) - C3*sum(XD,2) - 2*c/(a^2-c^2)*sw;
        zeros(n3,1)]; % (n1+n2+n3) x 1
    
    tmat = [1 1 1;
        1 1 0;
        0.5 1 0];
    
    f = -[Fa, Fb, Fc]*tmat/2/M; % (n1+n2+n3) x 3  dloglikelihood
end

if der == 2
    C1 = (a^2+b^2)/(a^2-b^2)^2; C2 = 2*a*b/(a^2-b^2)^2;
    C3 = (a^2+c^2)/(a^2-c^2)^2; C4 = 2*a*c/(a^2-c^2)^2;
    
    FaM = -C1*SM + C2*XM + 2*a/(a^2-b^2);
    Fb = C2*SM - C1*XM - 2*b/(a^2-b^2); % n1 x M
    
    FaD = -C3*SD + C4*XD + 2*a/(a^2-c^2);
    Fc = C4*SD - C3*XD - 2*c/(a^2-c^2); % n2 x M
    
    FaI = -1/a^2*SMD + 1/a;
    
    FaM2 = sum(FaM.^2); FaMb = sum(FaM.*Fb); Fb2 = sum(Fb.^2); % 1 x M
    fM = [mean(FaM2.*w),mean(FaMb.*w); mean(FaMb.*w),mean(Fb2.*w)];
    Mmat = [1 1;1 1;1 0]; 
    fM = Mmat*fM*Mmat';
    
    FaD2 = sum(FaD.^2); FaDc = sum(FaD.*Fc); Fc2 = sum(Fc.^2); % 1 x M
    fD = [mean(FaD2.*w),mean(FaDc.*w); mean(FaDc.*w),mean(Fc2.*w)];
    Dmat = [1 0.5;1 1;1 0]; 
    fD = Dmat*fD*Dmat';
    
    FaI2 = sum(FaI.^2); 
    fI = mean(FaI2.*w)*ones(3);
    
    f = (fM+fD+fI)/4; 
end

end


