function x = cg(A,b,setup)

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

xk = setup.x0;
rk1 = b-A*xk;
rktrk1 = rk1'*rk1;
pk = rk1;
for k = 1:setup.maxits
    Apk = A*pk;
    a = rktrk1/(pk'*Apk);
    xk = xk+a*pk;
    rk = rk1-a*Apk;
    rktrk = rk'*rk;
    pk = rk+(rktrk/rktrk1)*pk;
    if setup.showComments;fprintf('norm(b-A*x) = %e\n',rktrk);end
    if sqrt(rktrk) < setup.tol
        break;
    end
    rk1 = rk;
    rktrk1 = rktrk;
end

% Output a few things
if(rktrk < setup.tol)
    fprintf('CG CONVERGED at iteration %i. \nnorm(b-A*x) = %e\n',k,norm(b-A*xk));
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
%		precond:		 false
%		L:		 0
%		U:		 0
%		M:		 0
%		x0:		 0
function setup = ensureSetup(setup,b)
% This function ensures the correct setup of your structure

names = fieldnames(setup);
trueNames = {'maxits','tol','showComments','precond','L','U','M','x0'};
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

if ~isfield(setup, 'normCalc')
	setup.normCalc = round(setup.maxits/100);
end


end