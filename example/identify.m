% IDENTIFY.M

function [B0inv]=identify(A,P)

K=size(P,1);
A1=A(1:K,1:K); A2=A(1:K,K+1:2*K); A3=A(1:K,2*K+1:3*K); A4=A(1:K,3*K+1:4*K);

% Horizon 0: Structural impact multiplier matrix
L0=P;

% Horizon infinity: Long-run multiplier
Linf=inv(eye(K)-A1-A2-A3-A4)*P;       

% Generate matrix of IRFs for all relevant horizons
L=[L0;Linf];

% Zero restrictions for monetary policy shock and AD shock
z1=[1 0 0 0 0 0; 0 0 0 1 0 0];
z2=[0 0 0 1 0 0];

% Algorithm j=1,2,3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First iteration: j=1
Q1=[z1*L];

% Q1' = Q * T, where Q is orthogonal and T is upper triangular. 
[Q,R]=qr(Q1');

% Choose q1 to be the last row of Q. 
q1=Q(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second iteration
Q2=[z2*L; q1'];

% Q2' = Q * T, where Q is orthogonal and T is upper triangular. 
[Q,R]=qr(Q2');
q2=Q(:,end);

% Third iteration
Q3=[q1';q2'];
[Q,R]=qr(Q3');
q3=Q(:,end);

% Assemble structural impact multiplier matrix
B0inv=P*[q1 q2 q3];

% Flip signs of columns as appropriate to ensure positive diagonal elements in B0inv
for j=2:K
     if B0inv(j,j)<0
           B0inv(:,j)=-B0inv(:,j);
     end;
end;

% Normalize so a MP shock raises the Federal Funds rate
if B0inv(2,1)<0
    B0inv(:,1)=-B0inv(:,1);
end;    

% Normalize so a positive AS shock raises real GNP
if B0inv(1,3)<0
    B0inv(:,3)=-B0inv(:,3);
end;

