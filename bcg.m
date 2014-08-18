%% Block Conjugate Gradient
% A block conjugate gradient algorithm to solve multiple RHSs at once.
% Where A is a symmetric n*n PD matrix and B is a n*m matrix with the
% multiple RHSs in the columns.
%
% [X] = bcg(A,B);
%
% A setup structure can be passed in to switch between internal methods.
%
%  [X] = bcg(A,B,setup);
%
% The setup structure has fields with defaults:
%
% SETUP:
%		maxits:		 10
%		tol:		 1.000000e-06
%		rankTol:	 1.000000e-06
%		showComments:		 false
function [X k rnk nrm_rk] = bcg(A,B,setup)


% Set up the arguments
if nargin <2
    disp('You must include at least 2 arguments: A and b');
    X = ensureSetup(struct());
    return;
elseif nargin <3
    setup = ensureSetup(struct());
else
    setup = ensureSetup(setup);
end

% set up initial residual R0 = B;
[Rk1 r e] = qr(B,0);%With column pivoting
rnk = 0;
re = e;% reverse permutation vector for the end
for i = 1:length(e);
    if abs(r(i,i)) > setup.rankTol
        rnk = rnk+1;
    end
    re(e(i)) = i;
end
Rk1 = Rk1(:,1:rnk);%ensure economy sized
r = r(1:rnk,:);


Xk = zeros(size(Rk1));
Pk = Rk1;
p = 1:size(Rk1,2);% just an index
for k = 1:setup.maxits
    APk = A*Pk(:,p);
    a = (Pk(:,p)'*APk)\(Rk1(:,p)'*Rk1(:,p));
    Xk(:,p) = Xk(:,p)+Pk(:,p)*a;
    Rk(:,p) = Rk1(:,p)-APk*a;
    % The Removal of dependent columns using a column pivoting QR algorithm.
    % is not done here. But could be!
    % Instead:
    % Check the norm of each new residual.
    for n = length(p):-1:1
        nrm(p(n)) = norm(Rk(:,p(n)));
        if(nrm(p(n)) < setup.tol)
            % If the norm of one is below tolerance: remove it.
            p(n) = [];
        end
    end
    nrm_rk = max(nrm);
    if setup.showComments;fprintf('norm(b-A*x) = %e\n',nrm_rk);end
    if nrm_rk < setup.tol || isempty(p)
        break;
    end
    b = (Rk1(:,p)'*Rk1(:,p))\(Rk(:,p)'*Rk(:,p));
    Pk(:,p) = (Rk(:,p)+Pk(:,p)*b);
    Rk1(:,p) = Rk(:,p);
end

X = Xk*r;
X = X(:,re);
% Output a few things
nrm_rk = norm(B-A*X)/norm(B,'fro');% do one explicit calculation of the norm
if(nrm_rk < setup.tol)
    fprintf('CG CONVERGED at iteration %i. \nnorm(B-A*X)/norm(B,''fro'') = %e\n',k,nrm_rk);
else
    fprintf('\n\n*****************************\n\nminres DID NOT converge by iteration %i. \nnorm(B-A*X)/norm(B,''fro'') = %e\n\n*****************************\n\n',k,nrm_rk);
end


end

%% Ensure Setup Structure is Correct
% A setup structure can be passed in to switch between internal methods.
% The setup structure has fields:
%
% SETUP:
%		maxits:		 10
%		tol:		 1.000000e-06
%		rankTol:	 1.000000e-06
%		showComments:		 false
function setup = ensureSetup(setup)
% This function ensures the correct setup of your structure

names = fieldnames(setup);
trueNames = {'maxits','tol','rankTol','showComments'};
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
if ~isfield(setup, 'rankTol')
	setup.rankTol = 1.000000e-06;
end
if ~isfield(setup, 'showComments')
	setup.showComments = false;
end

end