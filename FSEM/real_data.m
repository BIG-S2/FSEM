clear; clc;

addpath('src')

load('data/demo.mat');

rng(1234567);
%% 1. Preprocessing

[MZtp1,~,DZtp1,~,MDZti] = FSEM_index(famid,zygosity); % famid and zygosity: N x 1 vector
n = sum(MZtp1) + sum(DZtp1) + sum(MDZti);

age_MRI = (age_MRI-mean(age_MRI))/std(age_MRI);
gender = (gender-mean(gender))/std(gender);
X = [ones(length(age_MRI),1),gender,age_MRI]; % X: N x 3 matrix - covariates matrix

FA = csvread('data/FA_GENU.csv', 1, 1)'; % FA: N x M matrix
FA(end,:) = [];
FA(:, any(isnan(FA), 1)) = [];

MD = csvread('data/MD_GENU.csv', 1, 1)'; % MD: N x M matrix
MD(end,:) = [];
MD(:, any(isnan(MD), 1)) = []; 

meanFA = mean(FA(:));
meanMD = mean(MD(:));
maxMean = max([meanFA, meanMD]);
FA = FA*(maxMean/meanFA);
MD = MD*(maxMean/meanMD);

sM = linspace(0, 1, size(FA,2));  % 1 x M vector
h = 0.02;
%% 2. FSEM for FA data
%% 2.1 Maximum likelihood estimation and local test
[FA_paramsA, FA_paramsN] = FSEM_mle(FA, X, famid, zygosity);
[FA_LRT_vec, FA_p_vec] = FSEM_lrt(FA_paramsA, FA_paramsN, FA, X, famid, zygosity);
%% 2.2 Weighted maximum likelihood estimation and local test
[FA_betaA, FA_log_WsigmaA, FA_log_WsigmaN] = FSEM_wmle(FA, X, famid, zygosity, sM, h);
FA_R = FA - X*FA_betaA;
[FA_WLRT_vec, FA_Wp_vec] = FSEM_wlrt(FA_log_WsigmaA, FA_log_WsigmaN, FA_R, famid, zygosity, sM, h);
%% 2.3 Global test
FA_p = FSEM_global_test(FA_R, famid, zygosity, 1e3);
%% 2.4 Covariance structure estimation
FA_ob = FSEM_cov_ob(FA_R, famid, zygosity, sM);
[FA_SigmaA_est, FA_SigmaC_est, FA_SigmaEG_est] = ...
    FSEM_cov(FA_R, famid, zygosity, sM, FA_ob);
%% 3. FSEM for MD data 
%% 3.1 Maximum likelihood estimation and local test
[MD_paramsA, MD_paramsN] = FSEM_mle(MD, X, famid, zygosity);
[MD_LRT_vec, MD_p_vec] = FSEM_lrt(MD_paramsA, MD_paramsN, MD, X, famid, zygosity);
%% 3.2 Weighted maximum likelihood estimation and local test
[MD_betaA, MD_log_WsigmaA, MD_log_WsigmaN] = FSEM_wmle(MD, X, famid, zygosity, sM, h);
MD_R = MD - X*MD_betaA;
[MD_WLRT_vec, MD_Wp_vec] = FSEM_wlrt(MD_log_WsigmaA, MD_log_WsigmaN, MD_R, famid, zygosity, sM, h);
%% 3.3 Global test
MD_p = FSEM_global_test(MD_R, famid, zygosity, 1e3);
%% 3.4 Covariance structure estimation
MD_ob = FSEM_cov_ob(MD_R, famid, zygosity, sM);
[MD_SigmaA_est, MD_SigmaC_est, MD_SigmaEG_est] = ...
    FSEM_cov(MD_R, famid, zygosity, sM, MD_ob);

%% 4. Confidence bands for beta
M = size(MD,2);
FA_paramsA_se = zeros(6,M); 
MD_paramsA_se = zeros(6,M);

for kk = 1:M
    tmp = hessian(@(params) FSEM_lla(params,FA(:,kk),X,famid,zygosity), ...
        FA_paramsA(:,kk));
    tmp = diag(inv(tmp)); 
    FA_paramsA_se(:,kk) = sqrt(2*tmp).*[exp(FA_paramsA(1:3,kk));ones(3,1)];
    
    tmp = hessian(@(params) FSEM_lla(params,MD(:,kk),X,famid,zygosity), ...
        MD_paramsA(:,kk));
    tmp = diag(inv(tmp)); 
    MD_paramsA_se(:,kk) = sqrt(2*tmp).*[exp(MD_paramsA(1:3,kk));ones(3,1)];
end

save real_data.mat

%% 5. Plot figures
eps = 1e-6;
%% 5.1 Local p values figure
figure; 
FA_params_LRT = exp(FA_paramsA(1:3,:));
FA_params_WLRT = exp(FA_log_WsigmaA(1:3,:));
FA_herit_LRT = FA_params_LRT(1,:)./sum(FA_params_LRT);
FA_herit_WLRT = FA_params_WLRT(1,:)./sum(FA_params_WLRT);

subplot(2,2,1); plot(FA_herit_LRT,'r'); hold on; 
plot(FA_herit_WLRT,'b'); hold on;
legend('LRT','WLRT')
title('FA: Estimated heritability along genu fiber tract','FontSize',10); 
xlabel('arclength-(a)','FontSize',10);
subplot(2,2,3); plot(-log10(max(FA_p_vec,eps)),'r'); hold on; 
plot(-log10(max(FA_Wp_vec,eps)),'b'); 
line([0,100],[-log10(0.05) -log10(0.05)],'Color','k');
legend('LRT','WLRT');
title('FA: -log10(p) values along genu fiber tract','FontSize',10); 
xlabel('arclength-(c)','FontSize',10);

MD_params_LRT = exp(MD_paramsA(1:3,:));
MD_params_WLRT = exp(MD_log_WsigmaA(1:3,:));
MD_herit_LRT = MD_params_LRT(1,:)./sum(MD_params_LRT);
MD_herit_WLRT = MD_params_WLRT(1,:)./sum(MD_params_WLRT);

subplot(2,2,2); plot(MD_herit_LRT,'r'); hold on; 
plot(MD_herit_WLRT,'b'); hold on; 
legend('LRT','WLRT')

title('MD: Estimated heritability along genu fiber tract','FontSize',10); 
xlabel('arclength-(b)','FontSize',10);
subplot(2,2,4); plot(-log10(max(eps,MD_p_vec)),'r'); hold on; 
plot(-log10(max(eps,MD_Wp_vec)),'b'); 
line([0,100],[-log10(0.05) -log10(0.05)],'Color','k');
legend('LRT','WLRT');
title('MD: -log10(p) values along genu fiber tract','FontSize',10); 
xlabel('arclength-(d)','FontSize',10);

%% 5.2 Varying coefficient figure
figure; 
subplot(2,3,1); plot(FA_paramsA(4,:),'b'); xlabel('arclength-(a)','FontSize',10);
ylabel('95% confidence band', 'FontSize', 10);
hold on; plot(FA_paramsA(4,:)+1.96*FA_paramsA_se(4,:),'--r'); 
         plot(FA_paramsA(4,:)-1.96*FA_paramsA_se(4,:),'--r');
         title('FA: intercept','FontSize',14);
subplot(2,3,2); plot(FA_paramsA(5,:)); xlabel('arclength-(b)','FontSize',10);
ylabel('95% confidence band', 'FontSize', 10);
hold on; plot(FA_paramsA(5,:)+1.96*FA_paramsA_se(5,:),'--r'); 
         plot(FA_paramsA(5,:)-1.96*FA_paramsA_se(5,:),'--r');
         title('FA: gender','FontSize',14);
subplot(2,3,3); plot(FA_paramsA(6,:)); xlabel('arclength-(c)','FontSize',10);
ylabel('95% confidence band', 'FontSize', 10);
hold on; plot(FA_paramsA(6,:)+1.96*FA_paramsA_se(6,:),'--r'); 
         plot(FA_paramsA(6,:)-1.96*FA_paramsA_se(6,:),'--r');
         title('FA: age','FontSize',14);
subplot(2,3,4); plot(MD_paramsA(4,:)); xlabel('arclength-(d)','FontSize',10);
ylabel('95% confidence band', 'FontSize', 10);
hold on; plot(MD_paramsA(4,:)+1.96*MD_paramsA_se(4,:),'--r'); 
         plot(MD_paramsA(4,:)-1.96*MD_paramsA_se(4,:),'--r');
         title('MD: intercept','FontSize',14);
subplot(2,3,5); plot(MD_paramsA(5,:)); xlabel('arclength-(e)','FontSize',10);
ylabel('95% confidence band', 'FontSize', 10);
hold on; plot(MD_paramsA(5,:)+1.96*MD_paramsA_se(5,:),'--r'); 
         plot(MD_paramsA(5,:)-1.96*MD_paramsA_se(5,:),'--r');
         title('MD: gender','FontSize',14);
subplot(2,3,6); plot(MD_paramsA(6,:)); xlabel('arclength-(f)','FontSize',10);
ylabel('95% confidence band', 'FontSize', 10);
hold on; plot(MD_paramsA(6,:)+1.96*MD_paramsA_se(6,:),'--r'); 
         plot(MD_paramsA(6,:)-1.96*MD_paramsA_se(6,:),'--r');
         title('MD: age','FontSize',14);

%% 5.3 Covariance structure figure
FA_lb = min([min(FA_SigmaA_est(:)), min(FA_SigmaC_est(:)), min(FA_SigmaEG_est(:))]); 
FA_ub = max([max(FA_SigmaA_est(:)), max(FA_SigmaC_est(:)), max(FA_SigmaEG_est(:))]); 
MD_lb = min([min(MD_SigmaA_est(:)), min(MD_SigmaC_est(:)), min(MD_SigmaEG_est(:))]); 
MD_ub = max([max(MD_SigmaA_est(:)), max(MD_SigmaC_est(:)), max(MD_SigmaEG_est(:))]); 

set(gcf,'units','points','position',[10,10,1200,500])
subplot(2,3,1); imagesc(FA_SigmaA_est, [FA_lb, FA_ub]); 
colorbar; title('FA: Estimated Genetic Covariance','FontSize',10);
xlabel('(a)','FontSize',10);
subplot(2,3,2); imagesc(FA_SigmaC_est, [FA_lb, FA_ub]); 
colorbar; title('FA: Estimated Common Environmental Covariance','FontSize',10);
xlabel('(b)','FontSize',10);
subplot(2,3,3); imagesc(FA_SigmaEG_est, [FA_lb, FA_ub]);
colorbar; title('FA: Estimated Unique Environmental Covariance','FontSize',10);
xlabel('(c)','FontSize',10);

subplot(2,3,4); imagesc(MD_SigmaA_est, [MD_lb, MD_ub]); 
colorbar; title('MD: Estimated Genetic Covariance','FontSize',10);
xlabel('(d)','FontSize',10);
subplot(2,3,5); imagesc(MD_SigmaC_est, [MD_lb, MD_ub]); 
colorbar; title('MD: Estimated Common Environmental Covariance','FontSize',10);
xlabel('(e)','FontSize',10);
subplot(2,3,6); imagesc(MD_SigmaEG_est, [MD_lb, MD_ub]);
colorbar; title('MD: Estimated Unique Environmental Covariance','FontSize',10);
xlabel('(f)','FontSize',10);