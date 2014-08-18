function [x] = MRes(A,b,setup)

% Set up the arguments
if nargin <2
    disp('You must include at least 2 arguments: A and b');
    x = ensureSetup(struct(),b);
    return;
elseif nargin <3
    setup = ensureSetup(struct(),b);
else
    setup = ensureSetup(setup,b);
end



% Set up Arguments
s = zeros(setup.maxits,1);c = s;%Givens rotations
% Allocate space for T
% There are only four entries per column.
% This has the main diag on 3, with one subdiag
% T is changed in place to R, and the last row of T can be ignored
% To look at T use:
%
%       T = spdiags(T(1:3,1:k)',[2 1 0],k,k);
%
T = zeros(4,setup.maxits);
% Set up Q and P
% You will only need 3 at a time, so cycle through them
Q = zeros(length(b),3);
P = Q;
% Set aside some space for the tridiagonal elements of T
alp = zeros(setup.maxits,1);% where alpha is the main diag
B = zeros(setup.maxits,1);%   and B is the first super- and sub-diagonal.


xk = setup.x0;
g = zeros(setup.maxits+1,1);%Unit Vector

if setup.precond%left preconditioning!
    %if you are using preconditioning, do a quick forward and back
    %solve to find your r_0* = M^-1 * r_0
    Q(:,2) = setup.U\(setup.L\(b-A*setup.x0));%
else
    Q(:,2) = (b-A*setup.x0);%(unnormalized)
end
% initiate g(1) to norm(r_0)
g(1) = norm(Q(:,2));
Q(:,2) = Q(:,2)./g(1);%Now normalize q_1

% Loop Until Convergence or Max Iterations
for k = 1:setup.maxits
    % If you have preconditioning
    if setup.precond
        z = setup.U\(setup.L\(A*Q(:,2)));
    else
        z = A*Q(:,2);
    end
    %Where Aqj =B(j?1)q(j?1) +alpha(j)q(j) +B(j)q(j+1)
    alp(k) = Q(:,2)'*z;%hit it with q_j, to get your alpha coefficient
    z = z - alp(k)*Q(:,2);%subtract off first two terms (above)
    if k > 1%first term is 0 if j=1 because B(0)=0 and q(:,0) = 0
        z = z - B(k-1)*Q(:,1);
    end
    B(k) = norm(z);%normalize to get B(j) because ||q(j+1)|| = 1
    if B(k) == 0;
        fprintf('Terminated at iteration %i, B(%i) == 0\n',it,k);
    else
        %if B(j) ~= 0 divide z by it to get your new orth. vec
        Q(:,3) = z./B(k);
    end
    
    T(3,k) = alp(k);
    if k > 1
        T([2,4],k) = [B(k-1); B(k)];
    else
        T(4,k) = B(k);
    end
    % Apply all the previous rotations to the k-th column of T.
    % this is a simplified matrix-vector product U'*t_k
    % Note that you only need to do the last 2 rotations.
    if k > 2
        T(1:2,k) = [c(k-2) s(k-2);-s(k-2) c(k-2)]*T(1:2,k);
    end
    if k > 1
        T(2:3,k) = [c(k-1) s(k-1);-s(k-1) c(k-1)]*T(2:3,k);
    end
    % Solve for the current Givens Rotation
    r = norm(T(3:4,k));
    s(k) = T(4,k)/r;
    c(k) = T(3,k)/r;
    % Apply the rotation to the sub-diagonal element in column k
    T(3,k) =r;% c(k).*T(3,k) + s(k).*T(4,k);
%   T(4,k) =   -s(k).*T(3,k) + c(k).*T(4,k); %This is zero.
    % Apply the rotation to g
    g(k+1) = -s(k).*g(k);
    g( k ) =  c(k).*g(k);
    
    
    % Update P = QR^-1
    P(:,3) = (Q(:,2)-T(2,k).*P(:,2)-T(1,k).*P(:,1))./T(3,k);
    
    xk = xk + g(k).*P(:,3);
    Q = Q(:,[2 3 1]);
    P = P(:,[2 3 1]);
    
    % Approximate the Norm
    % Approximate (due to rounding errors) the norm(r_k) by g(k+1)
    nrm_rk = abs(g(k+1));
    if setup.showComments;fprintf('norm(b-A*x) = %e\n',nrm_rk);end
    
    % Stopping Criteria
    if(nrm_rk < setup.tol)
        break;
    end
end %loop
% Output a few things
if(nrm_rk < setup.tol)
    fprintf('minres CONVERGED at iteration %i. \nnorm(b-A*x) = %e\n',k,norm(b-A*xk));
else
    fprintf('\n\n*****************************\n\nminres DID NOT converge by iteration %i. \nnorm(b-A*x) = %e\n\n*****************************\n\n',k,norm(b-A*xk));
end
x = xk;



end
%% Ensure Setup Structure is Correct
% A setup structure can be passed in to switch between internal methods.
% The setup structure has fields:
%
% SETUP:
%		maxits:		 10
%		tol:		 1.000000e-06
%		showComments:		 false
%		record:		 false
%		precond:		 false
%		L:		 0
%		U:		 0
%		M:		 0
%		x0:		 0
function setup = ensureSetup(setup,b)
% This function ensures the correct setup of your structure

names = fieldnames(setup);
trueNames = {'maxits','tol','showComments','record','precond','L','U','M','x0'};
for i =1:length(names)
	if ~any(strcmp(names{i},trueNames))
		warning('Setup:FieldNR','Field not recognized: ''%s''',names{i});
	end
end

if ~isfield(setup, 'maxits')
	setup.maxits = 10;
end

if ~isfield(setup, 'tol')
	setup.tol = 1.000000e-06;
end

if ~isfield(setup, 'showComments')
	setup.showComments = false;
end

if ~isfield(setup, 'record')
	setup.record = false;
end

if ~isfield(setup, 'precond')
	setup.precond = false;
end

if ~isfield(setup, 'L')
	setup.L = 0;
end

if ~isfield(setup, 'U')
	setup.U = 0;
end

if ~isfield(setup, 'M')
	setup.M = 0;
end

if ~isfield(setup, 'x0')
	setup.x0 = zeros(size(b));
end

end