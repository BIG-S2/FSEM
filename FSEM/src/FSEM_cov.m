function [SigmaA_est, SigmaC_est, SigmaE_est] = FSEM_cov(R, famid, zygosity, sM, h, null)
% FSEMcov computes the global variation matrix
% 
% argin
% R:        N x M matrix; residual matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector
% sM:       1 x M vector
% h:        scalar; global bandwidth
% null:     bool; estimation under null or not
% 
% argout:
% SigmaA_est: M x M matrix; global genetic variation
% SigmaC_est: M x M matrix; global environment variation
% SigmaE_est: M x M matrix; global error variation
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA

if nargin < 6
    null = false;
end

[N, M] = size(R);

SigmaA_est = zeros(M, M);
SigmaC_est = zeros(M, M);
SigmaE_est = zeros(M, M);

[MZtp1, MZtp2, DZtp1, DZtp2, ~] = FSEM_index(famid, zygosity);

S0 = R'*R/N;
S1 = (R(MZtp1,:)'*R(MZtp2,:) + R(MZtp2,:)'*R(MZtp1,:))/2/sum(MZtp1);
S2 = (R(DZtp1,:)'*R(DZtp2,:) + R(DZtp2,:)'*R(DZtp1,:))/2/sum(DZtp1);

for ss = 1:M
    for tt = ss:M
        vecs = sM-sM(ss);
        vect = sM-sM(tt);
        
        Wmat = 1/h^2*ker(vecs'/h)*ker(vect/h);
        W = sum(Wmat(:))-sum(diag(Wmat));
        
        Sw0mat = Wmat.*S0;
        Sw1mat = Wmat.*S1;
        Sw2mat = Wmat.*S2;
        
        Sw0 = sum(Sw0mat(:))-sum(diag(Sw0mat));
        Sw1 = sum(Sw1mat(:))-sum(diag(Sw1mat));
        Sw2 = sum(Sw2mat(:))-sum(diag(Sw2mat));
        if (null)
            SigmaC_est(ss, tt) = (Sw1+Sw2)/2/W;
            SigmaE_est(ss, tt) = (2*Sw0-Sw1-Sw2)/2/W;
        else
            SigmaA_est(ss,tt) = 2*(Sw1-Sw2)/W;
            SigmaC_est(ss,tt) = (-Sw1+2*Sw2)/W;
            SigmaE_est(ss,tt) = (Sw0-Sw1)/W;
        end
    end
end

SigmaA_est = (SigmaA_est+SigmaA_est'-diag(diag(SigmaA_est)));
SigmaC_est = (SigmaC_est+SigmaC_est'-diag(diag(SigmaC_est)));
SigmaE_est = (SigmaE_est+SigmaE_est'-diag(diag(SigmaE_est)));

end



