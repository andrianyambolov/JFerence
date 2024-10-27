function lnf = irfpdf2cum(A0,Aplus,C,H,n,nutilde,Omegatilde,p,Phitilde,Psitilde)

% Purpose:
% This code computes the log of the joint posterior density of structural impulse responses up to constant 
% input : 
% A0,Aplus   : nxn and (np+1)xn matrices of structural VAR parameters in Arias et al. (2018).
% C          : (n(n-1)/2)x(2n^2) matrix of constants for short-run and long-run identification restrictions
% H          : Maximum horizon
% n          : number of variables
% nutilde    : Degrees of freedom for the inverse-Wishart distribution. nuT=nu0+T in Uhlig (2005).
% Omegatilde : (np+1)x(np+1) matrix.  inv(NT)=inv(N0+X'*X) in Uhlig (2005)  
% p          : number of lags
% Phitilde   : nxn matrix. Y'*Y+Phibar+Psibar'*inv(Omegabar)*Psibar-Psitilde'*inv(Omegatilde)*Psitilde in Arias et al. (2018).  
%              nuT*ST in Uhlig (2005).
% Psitilde   : (np+1)xn matrix.  BbarT in Uhlig (2005) with the last element of Xt being one (i.e., the last row of Psitilde is the intercept term.) 
% output:
% lnf     : Log of the joint posterior density up to constant

Kn=[]; % Commutation matrix
for i=1:n
    for j=1:n
        E      = zeros(n,n);
        E(i,j) = 1;
        Kn     = [Kn reshape(E,n*n,1)];
    end
end;
Dn = []; % Duplication matrix
mask = tril(true(n,n));
for i=1:n
    for j=1:n
        E      = zeros(n,n);
        E(i,j) = 1;
        E(j,i) = 1;
        Dn = [Dn;E(mask)']; 
    end
end
Dnplus = inv(Dn'*Dn)*Dn';
L = min(p,H);
 
% log of NIW
invA0      = inv(A0);
Sigma      = invA0'*invA0; %A0 is defined as in Arias et al. (2018)
B          = Aplus*invA0;
invB1      = eye(n);
for i=1:p
    invB1 = invB1-B(n*(i-1)+1:n*i,:);
end
invB1 = inv(invB1); 
Omegatilde = Omegatilde(1:n*p,1:n*p);
lnf = -0.5*n*p*log(abs(det(Sigma)))+0.5*n*log(abs(det(Omegatilde)))-0.5*(reshape(B(1:n*p,:)-Psitilde(1:n*p,:),n*n*p,1))'*kron(inv(Sigma),Omegatilde)*reshape(B(1:n*p,:)-Psitilde(1:n*p,:),n*n*p,1);
lnf = lnf-0.5*(nutilde-n-1)*log(abs(det(Sigma)))+0.5*nutilde*log(abs(det(Phitilde)))-0.5*trace(Phitilde*inv(Sigma));

% Compute reduced-form impulse responses
M = [eye(n) zeros(n,n*(p-1))]';
F = [B(1:n*p,:)';eye(n*(p-1)) zeros(n*(p-1),n)]; % Companion matrix 
Phi = eye(n);
for i = 1:H
    Phi = [Phi M'*(F^i)*M];
end
s = zeros(n,1);
s(1,1) = 1;
S = diag(s); 
J1 = eye((H+1)*n*n);
for i=2:H
    for j=1:i-1
        J1((i-1)*n*n+1:i*n*n,(j-1)*n*n+1:j*n*n) = kron(eye(n),S);
    end
end
J2J3 = zeros((H+1)*n*n,n*n*(L+1));
J2J3(1:n*n,1:n*n) = -kron(invA0,invA0')*Kn;
for i=1:H
    J2J3(i*n*n+1:(i+1)*n*n,1:n*n) = -kron(invA0,Phi(:,i*n+1:(i+1)*n)*invA0')*Kn;
    if i<=L
        J2J3(i*n*n+1:(i+1)*n*n,i*n*n+1:(i+1)*n*n) = kron(invA0,eye(n)); 
    end
    if i>1
        invA0kronIX = zeros(n*n,n*n*p);
        for k=0:(i-1)
            invA0kronIX = invA0kronIX+kron(invA0*M'*(F')^k,Phi(:,(i-1-k)*n+1:(i-k)*n)');
        end
        invA0kronIX = invA0kronIX*kron(eye(p),Kn); 
        if i<=L
            J2J3(i*n*n+1:(i+1)*n*n,n*n+1:i*n*n) = invA0kronIX(:,1:(i-1)*n*n);
        else
            J2J3(i*n*n+1:(i+1)*n*n,n*n+1:end) = invA0kronIX;
        end
    end
end
D = zeros(n*n*(L+1),n*n*(L+1));
D(1:n*(n+1)/2,1:n*n) = Dnplus*(kron(invA0',Sigma)+kron(Sigma,invA0')*Kn);
N = zeros(n*n*(L+1),n*(n+1)/2+n*n*L);
CA  = C*[kron(invA0,invA0')*Kn; kron(invA0,invB1*invA0')*Kn];
CB  = zeros(n*(n-1)/2,n*n*L);
for i=1:L
    CB(:,1+n*n*(i-1):n*n*i) = C*[zeros(n*n,n*n);kron(invA0*invB1',invB1)];
end;
D(n*(n+1)/2+1:n*n,1:n*n) = CA;
D(n*n+1:n*n*(1+L),n*n+1:end) = eye(n*n*L);
N(1:n*(n+1)/2,1:n*(n+1)/2) = eye(n*(n+1)/2);
N(n*(n+1)/2+1:n*n,n*(n+1)/2+1:end) = CB;
N(n*n+1:n*n*(L+1),n*(n+1)/2+1:end) = eye(n*n*L);

J4 = inv(D)*N;
lnf = lnf+0.5*logdet(J4'*J2J3'*J1'*J1*J2J3*J4);
