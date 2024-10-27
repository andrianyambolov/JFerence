% IRFVAR.M

function [IRF]=irfvar(A,B0inv,p,h)

q=size(B0inv,1);
J=[eye(q,q) zeros(q,q*(p-1))];
IRF=reshape(J*A^0*J'*B0inv,q^2,1);

for i=1:h-1
	IRF=([IRF reshape(J*A^i*J'*B0inv,q^2,1)]);
end;

