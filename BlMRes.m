% A block minimum residual Code.
% With householder reflections and optimal storage.
%
% Uses a Ruhe variation to calculate lanczos.
function [X] = BlMRes(A,B,setup)

% Set up the arguments
if nargin <2
    disp('You must include at least 2 arguments: A and b');
    X = ensureSetup(struct(),B);
    return;
elseif nargin <3
    setup = ensureSetup(struct(),B);
else
    setup = ensureSetup(setup,B);
end



% Set up Arguments
p = rank(full(B));
w = zeros(p+1,setup.maxits);%This p,k vector holds all the reflection info
% Allocate space for T
% There are only four entries per column.
% This has the main diag on 3, with one subdiag
% T is changed in place to R, and the last row of T can be ignored
% To look at T use:
%
%       T = spdiags(T(1:3,1:k)',[2 1 0],k,k);
%
T = zeros(p*3+1,setup.maxits+p);
% p*2+1 is the diagonal
d = p*2+1;
symInd = sub2ind([size(T,1),inf],d-(1:p),1+(1:p));
symIndUpd = p*3+1;
R = zeros(p*3+1,setup.maxits);
% Set up Q and P
% You will only need 3 at a time, so cycle through them
Q = zeros(length(B),1+2*p);
P = Q;
cShift = [2:(2*p+1) 1];

log_nrm_rk = inf;
logTol = log10(setup.tol);
cnt = 1;


Xk = setup.x0;
G = zeros(setup.maxits+1,size(B,2));

if setup.precond%left preconditioning!
    %if you are using preconditioning, do a quick forward and back
    %solve to find your R_0* = M^-1 * R_0
    [q r] = qr(setup.U\(setup.L\(B-A*X0)),0);
else
    [q r] = qr(B-A*setup.x0,0);
end

Q(:,p+1:end-1) = q(:,1:p);
G = G.*0;
G(1:p,1:size(B,2)) =  r(1:p,1:size(B,2));
skip = 0;
% Loop Until Convergence or Max Iterations
for k = p:setup.maxits
    
    j = k-p+1;
    % If you have preconditioning
    if setup.precond
        z = setup.U\(setup.L\(A*Q(:,p+1)));
    else
        z = A*Q(:,p+1);
    end
    % Subtract off the info we know
    z = z - Q(:,1:p)*T((1:p)+p,j);
    % Create the info We don't know
    for i = p+1:(2*p)
        T(i+p,j) = Q(:,i)'*z;
        z = z - T(i+p,j).*Q(:,i);
    end
    T(i+p+1,j) = norm(z);
    T(symInd) = T(d+1:end,j);
    symInd = symInd + symIndUpd;
    if abs(T(i+p+1,j)) < setup.tol;
        fprintf('Rank Deficient\n');
        skip = skip + 1;
    end
    %if B(j) ~= 0 divide z by it to get your new orth. vec
    Q(:,2*p+1) = z./T(i+p+1,j);
    
    % QR Factorize H using Housholder Reflections
    %
    %           Qj = I - 2*w*w';
    %
    % Where:
    %           w = z./norm(z);
    %
    %           z = X(:,j); z(1) = z(1) + beta;
    %
    %           beta = sign(X(j,j))*norm(X(:,j));
    %
    R(:,j) = T(:,j);
    % Apply all the previous reflections to the k-th column of T.
    % this is a simplified matrix-vector product Qj'*R_k
    for ref = max(1,j-2*(p)):j-1
        rind = max(1,d-j+ref):d-j+ref+p;
        wind = length(rind)-1;
%         pad = zeros(p+1-length(rind),1);
%         v = 2*([pad;R(rind,j)]'*w(:,ref));%This is a scalar!
%         R(rind,j) = R(rind,j) - w(:,ref).*v;
        v = 2*(R(rind,j)'*w(end-wind:end,ref));%This is a scalar!
        R(rind,j) = R(rind,j) - w(end-wind:end,ref).*v;
    end

    % Solve for the current Householder Reflection in short form
    % We are getting wj in:   
    %
    %           Qj = I - 2*wj*wj';
    %
    
    beta = sign(R(d,j))*norm(R(d:end,j));% Go to the p+1 below the diag
    w(1,j) = beta+R(d,j);%z_j,j         Add beta to first element
    w(2:p+1,j) = R(d+1:d+p,j);%         Here there are only p elements
    w(:,j) = w(:,j)./norm(w(:,j));%     Normalize by the length

    % Apply the reflection to the sub-diagonal elements in column j
    % Using:
    % 
    %       (I - 2ww')X = X-wv'
    %
    %                 v = 2X'w
    %
    v = 2*R(d:end,j)'*w(:,j);%This is a scalar!
    R(d:end,j) = R(d:end,j) - w(:,j).*v;

    % Apply the reflection to G
    v = 2*G(j:j+p,:)'*w(:,j);%This is not a scalar.
    G(j:j+p,:) = G(j:j+p,:) - w(:,j)*v';
    
    
    % Update P = QR^-1
    P(:,2*p+1) = Q(:,p+1);
    % Do the forward solve on the transpose of P
    %    R'P' = Q'
    for f = 1:2*p
        P(:,2*p+1) = P(:,2*p+1) - R(f,j).*P(:,f);
    end
    P(:,2*p+1) = P(:,2*p+1)./R(d,j);
    
    if rem(j,p) == 0%Update all of X at one time.
        Xk = Xk + P(:,p+2:2*p+1)*G(j-p+1:j,:);
    end
    
    Q = Q(:,cShift);
    P = P(:,cShift);
    
    % Approximate the Norm
    % Only do this if you are close or evey p iterations
    if log_nrm_rk-logTol< 0.1 || rem(cnt-1,p) == 0
        nrm_rk = norm(G(j+1:end,:),'fro');
        log_nrm_rk = log10(nrm_rk);
    end
    cnt = cnt+1;
    if setup.showComments;fprintf('norm(b-A*x) = %e\n',nrm_rk);end
    
    % Stopping Criteria
    if(nrm_rk < setup.tol)
        break;
    end
    
end %loop
% Output a few things
if(nrm_rk < setup.tol)
    fprintf('minres CONVERGED at iteration %i. \nnorm(b-A*x) = %e\n',k,norm(B-A*Xk));
else
    fprintf('\n\n*****************************\n\nminres DID NOT converge by iteration %i. \nnorm(b-A*x) = %e\n\n*****************************\n\n',k,norm(B-A*Xk));
end
X = Xk;



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