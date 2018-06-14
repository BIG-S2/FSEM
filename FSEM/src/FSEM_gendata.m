function [y, X, famid, zygosity, eta] = FSEM_gendata(beta, SigmaA, SigmaC, SigmaEG, SigmaEL, nn)

n1 = nn(1); n2 = nn(2); n3 = nn(3);

N = 2*(n1 + n2) + n3; p = size(beta, 1); M = size(SigmaA, 1);

if length(SigmaEL) == 1
    SigmaEL = ones(M, 1)*SigmaEL;
end

famid = [reshape(repmat(1:(n1+n2),2,1),2*(n1+n2),1);
         reshape(n1+n2+(1:n3),n3,1)];

zygosity = [zeros(2*n1,1); ones(2*n2+n3,1)];

% correlation matrix for x2,...,xp, compound sysmetric with rho = 1/sqrt(2)
SigmaX = ones(p-1,p-1)/sqrt(2);
SigmaX(linspace(1,(p-1)^2,p-1)) = 1; % set diagonal to be 1

X = [ones(N,1), mvnrnd(zeros(1,p-1), SigmaX, N)];
% standardize predictors
X(:,2:end) = bsxfun(@minus, X(:,2:end), mean(X(:,2:end)));
X(:,2:end) = bsxfun(@rdivide, X(:,2:end), std(X(:,2:end)));

AA = mvnrnd(zeros(1,M), SigmaA, n1+3*n2+n3);
CC = mvnrnd(zeros(1,M), SigmaC, n1+n2+n3);
EG = mvnrnd(zeros(1,M), SigmaEG, N);

AG = [reshape(repmat(AA(1:n1,:)',2,1),M,2*n1)'; % MZ
    sqrt(0.5)*reshape(repmat(AA(n1+(1:n2),:)',2,1),M,2*n2)' + ...
    sqrt(0.5)*AA(n1+n2+(1:2*n2),:); % DZ
    AA(n1+3*n2+(1:n3),:)]; % IN

CE = [reshape(repmat(CC(1:(n1+n2),:)',2,1),M,2*(n1+n2))';
    CC(n1+n2+(1:n3),:)];

eta = real(AG + CE + EG);

y = real(X*beta + eta + bsxfun(@times, randn(N, M), sqrt(SigmaEL')));

end




