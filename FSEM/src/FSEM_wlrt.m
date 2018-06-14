function [WLRT_vec, Wp_vec] = ...
    FSEM_wlrt(log_WsigmaA, log_WsigmaN, R, famid, zygosity, sM, h, nb, method)
% FSEM_wlrt does pointwise test given parameter estimates via WLRT
%
% argin
% log_WsigmaA: 3 x M matrix
% log_WsigmaN: 2 x M matrix
% R:           N x M matrix
% famid:       N x 1 vector
% zygosity:    N x 1 vector
% sM:          1 x M vector; tract points
% h:           scalar or 1 x M vector
% nb:          number of bootstrap; 500 for simulation, 1e6 for real data
% method:      'approx' or 'exact', 
%              'exact' for real data, 'approx' for simulation data
%
% argout
% WLRT_vec: 1 x M vector; test statistics
% Wp_vec:   1 x M vector; p values
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA

if nargin < 9
    method = 'exact';
end

if nargin < 8
    method = 'exact';
    nb = 1e6;
end

loc = sM;

if nargin < 7
    h = 5*range(sM)/length(sM)*ones(1,length(loc));
end

K = length(loc);

if length(h) == 1
    h = h*ones(1, K);
end

[MZtp1, ~, DZtp1, ~, MDti] = FSEM_index(famid, zygosity);
n1 = sum(MZtp1); n2 = sum(DZtp1); n3 = sum(MDti);
n = n1 + n2 + n3;

ksi = randn(nb, n);

WLRT_vec = zeros(1, K);
Wp_vec = zeros(1, K);

for kk = 1:K
    % score matrix n x 3 
    J_n = FSEM_wlla([-100;log_WsigmaN(:,kk)], R, famid, zygosity, sM, ...
        loc(kk), h(kk), 1);
    J_n = J_n/sqrt(n); 
    % information matrix
    I_n = FSEM_wlla([-100;log_WsigmaN(:,kk)], R, famid, zygosity, sM, ...
        loc(kk), h(kk), 2)/n;
    [V, D] = eig(I_n);
    D = D(1,1);
    
    if strcmp(method, 'exact')
        % for real data sets        
        fn = FSEM_wlln(log_WsigmaN(:,kk), R, famid, zygosity, sM, ...
            loc(kk), h(kk));
        fa = FSEM_wlla(log_WsigmaA(:,kk), R, famid, zygosity, sM, ...
            loc(kk), h(kk));
        WLRT_vec(kk) = 2*max(0, fn-fa);
    else
        % for simulation data sets (this is asymptotic approximation of the above version)
        CC = sum(J_n, 1)*V;
        CC = CC(1);
        WLRT_vec(kk) = max(0, CC^2*(CC > 0)/D);
    end
    
    J_n = ksi*J_n; % nb x 3 matrix
    C1 = J_n*V;    % nb x 3 matrix
    C1 = C1(:,1);  % nb x 1 vector
    
    T_BB = max(0, C1.^2.*(C1*D > 0)/D);
    
    Wp_vec(kk) = mean(T_BB > WLRT_vec(kk));
end

end