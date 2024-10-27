% NIW_INPUTS.M

function [nutilde,Omegatilde,Phitilde,Psitilde]=niw_inputs(y, p)

% Conditional MLE/LS estimator for model with intercept
[t,n]=size(y); T=t-p;

Y=y(p:t,:);	
for i=1:p-1
 	Y=[Y y(p-i:t-i,:)];
end;
X=[Y(1:t-p,:) ones(T,1)];
Bhat=inv(X'*X)*X'*y(p+1:t,:);
B=Bhat';
Sigmahat = (y(p+1:t,:)-X*Bhat)'*(y(p+1:t,:)-X*Bhat)/T;
Ydep=y(p+1:t,:);

% Our NIW prior parameters
Bbar0=zeros(p*n+1,n);
nu0=0;
N0=zeros(n*p+1,n*p+1);
S0=zeros(n,n); 

% Our NIW posterior parameters
nuT   = T+nu0;
NT    = N0+X'*X;
BbarT = inv(NT)*(N0*Bbar0+X'*X*Bhat);
ST    = (nu0/nuT)*S0+(T/nuT)*Sigmahat+(1/nuT)*(Bhat-Bbar0)'*N0*inv(NT)*X'*X*(Bhat-Bbar0);

% Map into Arias et al. (2018) notation for NIW distribution
% Prior parameters
nubar=nu0;
Psibar=Bbar0;
%Omegabar=inv(N0); % Not well defined for our prior
Phibar=S0;

% Posterior parameters
nutilde=nuT;
Psitilde=BbarT;
Omegatilde=inv(NT);
Phitilde=Ydep'*Ydep + Phibar + Psibar'*N0*Psibar - Psitilde'*NT*Psitilde; % Note: Replaced inv(Omegabar) by N0 and inv(Omegatilde) by NT