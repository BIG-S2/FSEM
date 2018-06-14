function [ob] = FSEM_cov_ob(R, famid, zygosity, sM)
% FSEM_cov_ob computes optimal bandwidth via cross validation
%   for covariance matrix estimation
%
% argin
% R:        N x M matrix; residual matrix
% famid:    N x 1 vector
% zygosity: N x 1 vector
% sM:       1 x M vector
%
% argout
% ob:       optimal bandwidth
%
% Reference
%   S Luo, R Song, M Styner, JH Gilmore & H Zhu FSEM: Functional Structural
%         Equation Models for Twin Functional Data, JASA


    function [cvsse] = FSEM_cvsse(R, famid, zygosity, Sigma_aP_est, Sigma_cP_est, Sigma_eP_est)
        % FSEMcvsse computes the sum of squared errors via cross validation
        %   for given twin pairs or twin individual
        %
        % argin
        % R:            n x M matrix; residual matrix
        % famid:        n x 1 vector assume <sorted> from small to large 
        % zygosity:     n x 1 vector
        % Sigma_aP_est: M x M matrix; global genetic variation
        % Sigma_cP_est: M x M matrix; global environment variation
        %
        % argout
        % cvsse:        cross validation sum of squared errors
        
        cvsse = 0;
        [MZtp1, ~, DZtp1, ~, MDZti] = FSEM_index(famid, zygosity);
        ind = 1:size(R,1); indMZtp1 = ind(MZtp1); 
        indDZtp1 = ind(DZtp1); indMDZti = ind(MDZti);
        if (~isempty(indMZtp1))
            for k = 1:length(indMZtp1)
                tempA = ((R(indMZtp1(k),:)'*R(indMZtp1(k)+1,:) + ...
                    R(indMZtp1(k)+1,:)'*R(indMZtp1(k),:))/2- Sigma_aP_est - Sigma_cP_est - Sigma_eP_est).^2;
                cvsse = cvsse + sum(tempA(:)) - sum(diag(tempA));
            end
        end
        if (~isempty(indDZtp1))
            for k = 1:length(indDZtp1)
                tempA = ((R(indDZtp1(k),:)'*R(indDZtp1(k)+1,:) + ...
                    R(indDZtp1(k)+1,:)'*R(indDZtp1(k),:))/2- Sigma_aP_est - Sigma_cP_est - Sigma_eP_est).^2;
                cvsse = cvsse + sum(tempA(:)) - sum(diag(tempA));
            end
        end
        if (~isempty(indMDZti))
            for k = 1:length(indMDZti)
                tempA = (R(indMDZti(k),:)'*R(indMDZti(k),:)-Sigma_aP_est-Sigma_cP_est - Sigma_eP_est).^2;
                cvsse = cvsse + sum(tempA(:)) - sum(diag(tempA));
            end
        end
    end

n = famid(end);
nIndex = (1:n);
M = length(sM);
K = 5;% 5 fold cross validation
CVObj = cvpartition(n, 'Kfold', K);

% crude search
nh = 5;
xrange = range(sM);
hmin = xrange/M;
hmax = xrange*10/M;
vh = linspace(hmin, hmax, nh);
cvs = zeros(nh, 1);

for hh = 1:nh
    for kk = 1:K
        trainingIndex = ismember(famid, nIndex(CVObj.training(kk)));
        testIndex = ~trainingIndex;
        [Sigma_aP_est, Sigma_cP_est, Sigma_eP_est] = ...
            FSEM_cov(R(trainingIndex,:), famid(trainingIndex), zygosity(trainingIndex), sM, vh(hh));
        cvs(hh) = cvs(hh) + ...
            FSEM_cvsse(R(testIndex,:), famid(testIndex), zygosity(testIndex), Sigma_aP_est, Sigma_cP_est, Sigma_eP_est);
    end
end

[~, cvh] = min(cvs);
ob = vh(cvh);

% fine search
hmin = max(hmin, ob-0.02);
hmax = min(hmax, ob+0.02);
vh = linspace(hmin, hmax, nh);
cvs = zeros(nh, 1);

for hh = 1:nh
    for kk = 1:K
        trainingIndex = ismember(famid, nIndex(CVObj.training(kk)));
        testIndex = ~trainingIndex;
        [Sigma_aP_est, Sigma_cP_est, Sigma_eP_est] = ...
            FSEM_cov(R(trainingIndex,:), famid(trainingIndex), zygosity(trainingIndex), sM, vh(hh));
        cvs(hh) = cvs(hh) + ...
            FSEM_cvsse(R(testIndex,:), famid(testIndex), zygosity(testIndex), Sigma_aP_est, Sigma_cP_est, Sigma_eP_est);
    end
end

[~, cvh] = min(cvs);
ob = vh(cvh);

end