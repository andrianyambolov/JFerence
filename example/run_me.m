% Based on code by Lutz Kilian for the following paper:
% Kilian, L. and Inoue, A. (2022) Joint Bayesian Inference about Impulse 
% Responses in VAR Models. Journal of Econometrics, 231: 457-476.
% Adjusted by Andrian Yambolov, 13.04.2024

clc
clear
close all

rng(111) % random number generator seed
addpath(genpath(fullfile('..'))) % functions path

%% PARAMETERS

smpl_start = 'Q1-1939'; % sample start
smpl_end   = 'Q4-2015'; % sample end
p = 4; % lags
H = 36; % irf horizon
M=20000; % draws

%% DATA

raw = readtimetable('dataset.xlsx', 'Sheet', 'us');
est_time = datetime(smpl_start):calmonths(3):datetime(smpl_end);
data = raw(est_time,["log_rgov_exp_cap", "log_rgdp_cap", "log_rgov_rev_cap"]);
y = data{:,:};
[t,n] = size(y); % number observations and variables
T=t-p;

%% PRIORS

% reduced form NIW prior and posterior
Y=y(p:t,:);	
for i=1:p-1
 	Y=[Y y(p-i:t-i,:)];
end
X=[ones(T,1) Y(1:t-p,:)];
Bhat=inv(X'*X)*X'*y(p+1:t,:);
B=Bhat';
Sigmahat = (y(p+1:t,:)-X*Bhat)'*(y(p+1:t,:)-X*Bhat)/T;
Ydep=y(p+1:t,:);

% prior parameters
Bbar0=zeros(p*n+1,n);
nu0=0;
N0=zeros(n*p+1,n*p+1);
S0=zeros(n,n); 

%% ESTIMATION

disp('Start estimation ...')

% posterior parameters
nuT   = T+nu0;
NT    = N0+X'*X;
BbarT = inv(NT)*(N0*Bbar0+X'*X*Bhat);
ST    = (nu0/nuT)*S0+(T/nuT)*Sigmahat+(1/nuT)*(Bhat-Bbar0)'*N0*inv(NT)*X'*X*(Bhat-Bbar0);
EvecB = reshape(BbarT,n*(n*p+1),1); % Posterior mean of vec(B)

Sigmamat = [];
for i=1:M
    RANTR=chol(inv(ST))'*randn(n,nuT)/sqrt(nuT);
    Sigma=inv(RANTR*RANTR');
    Sigmamat = [Sigmamat; reshape(Sigma,1,n*n)];    
end

% Computing NIW inputs for density calculation 
[nutilde,Omegatilde,Phitilde,Psitilde]=niw_inputs(y, p);

irf = NaN(H, 3, M); % irf storage
% Evaluations of the posterior pdf at posterior draws
for j=1:M
    if ~mod(j, 5000)
        disp(['   - draw ', num2str(j) , '/', num2str(M)])
    end
    Sigma  = reshape(Sigmamat(j,1:n*n),n,n);
    VvecB = kron(Sigma,inv(NT));
    vecB   = EvecB+chol(VvecB)'*randn(n*(n*p+1),1);
    A      = (reshape(vecB,1+n*p,n))';
    B0inv      = chol(Sigma)';

    % Compute and store implied impulse response
    irf_draw = irfvar([A(:,2:p*n+1); eye((p-1)*n,(p-1)*n) zeros((p-1)*n,n)],B0inv,p,H);

    irf(:,:,j) = irf_draw(1:3,:)';

end


% save results
save('estimation.mat', 'data', 'irf');

