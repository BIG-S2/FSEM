function [pval, ts] = FSEM_global_test(R, famid, zygosity, nb)
% FSEM_global does global test FSEM: Functional Structural Equation Models
% for Twin Functional Data
% 
% argin
% R:        N x M matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector
% 
% argout
% pval: scalar P value
% ts:   scalar test statistic
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA


if nargin < 4
    nb = 1e3;
end

[MZtp1, MZtp2, DZtp1, DZtp2, ~] = FSEM_index(famid, zygosity);

% only randomize among twin pairs!!!
twins = MZtp1|MZtp2|DZtp1|DZtp2;
R = R(twins,:);
famid = famid(twins);
zygosity = zygosity(twins);
[MZtp1, MZtp2, DZtp1, DZtp2, ~] = FSEM_index(famid, zygosity);

N = size(R, 1);
M = size(R, 2);

sM = linspace(0,1,M);
knots = [min(sM), quantile_knots(sM,10), max(sM)];

bM = bsplineMat(sM, knots, 4);
pbM = (bM'*bM)\bM';

cM = R*pbM';
MZ = cM(MZtp1,:) - cM(MZtp2,:);
DZ = cM(DZtp1,:) - cM(DZtp2,:);
DZ = DZ'*DZ;
MZ = MZ'*MZ;
D = eig(DZ, MZ);
ts = max(D);

tss = zeros(nb,1);

for ii = 1:nb
    cM = R(randperm(N),:)*pbM';
    MZ = cM(MZtp1,:) - cM(MZtp2,:);
    DZ = cM(DZtp1,:) - cM(DZtp2,:);
    DZ = DZ'*DZ;
    MZ = MZ'*MZ;
    D = eig(DZ, MZ);
    tss(ii) = max(D);
end

pval = mean(tss > ts);

end