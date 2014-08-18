function [Q,H,R] = blockArnoldi(A,B,m,X0)

Q = zeros(length(B),m+1);
H = zeros(m+1,m);
[q R] = qr(B-A*X0,0);
p = size(q,2);
Q(:,1:p) = q;


for k = p:p+m
    j = k-p+1;
    z = A*Q(:,j);
    for i = 1:k
        H(i,j) = Q(:,i)'*z;
        z = z - H(i,j).*Q(:,i);
    end
    H(k+1,j) = norm(z);
    if H(k+1,j) == 0;
        fprintf('Terminated at k = %i, H(%i,%i) == 0\n',m,k+1,k);
        break;
    end
    Q(:,k+1) = z./H(k+1,j);
end
