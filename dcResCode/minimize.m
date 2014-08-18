function [x, k] = minimize( para, x, optTol, maxit, setup)


fnc = @para.getMisfit;
para.mref = x;

% Parse the options structure if it exists
if ~exist('setup','var')
    setup = ensureSetup(struct());
else
    setup = ensureSetup(setup);
end


% initialize BFGS vectors
Y       = zeros(length(x(:)) , setup.nbfgs);
S       = Y;
bfgscnt = 0;


H0 = para.GRADw'*para.GRADw;
% MM = @(z) tril(H)\(diag(H).*(triu(H)\z));
% Hinv = @(z) pcg(H,z,1e-6,200,MM);


% Evaluate the Function
    [f, g] = fnc( x );
% Start the minimization
for k = 1:maxit
    gold = g;
    %     % Do a BFGS update if we are not in the first iteration
    %         if k > 1 % not for the first one
    %             s = t*p;
    %             y = g_new-g;
    %             % Update the inverse of the hessian (Messy)
    %             % Use LM-BFGS
    %             H = H + (s'*y+y'*H*y)*(s*s')./((s'*y)^2)-(H*y*s'+s*(y'*H))/(s'*y);
    %
    %         end
    fprintf('%i)\t%e\t%e\t',k,norm(f),norm(g))
    p = bfgs(bfgscnt,S,Y,H0,-g);
%     p = Hinv(-g);

    % Scale p if needed.
    if norm(p,'inf') > setup.maxStepSize, 
        p = p/norm(p,'inf')*setup.maxStepSize; 
    end
    
    
    if norm(g) <= optTol;break;end
    
    [t numBacktrack] = backtrack(fnc,x,p,f,g,setup);%backtracking line search
    if numBacktrack >= 2
        setup.maxStepSize = setup.maxStepSize*0.75;
    elseif numBacktrack == 0
        setup.maxStepSize = min(1,setup.maxStepSize*2);
    end
    fprintf('\tmSS: %e\t',setup.maxStepSize);
        
    fprintf('\n');
    %update x
    x = x +t*p;
    if max(abs(t*p)) < setup.minStepSize
        disp('Converged based on size of model update.')
        break
    end
    if numBacktrack > setup.maxLSits
        disp('Not getting anywhere.')
        break
    end
    if setup.plotIt;montageArray(reshape(x,para.dims));
        caxis([-10 0]);drawnow;
    end
    
    [f, g] = fnc( x );
    
    % Update BFGS vectors
    yy = g - gold;
    ss = t*p;
    if yy'*ss > 0
        ktop      = mod(bfgscnt,setup.nbfgs)+1;
        Y(:,ktop) = yy;
        S(:,ktop) = ss;
        bfgscnt   = bfgscnt+1;
    else
        disp('   yTs < 0, skip update');
    end;
    
end

end



%% ARMIJO BACKTRACKING LINE SEARCH
%
% The backtracking algorithm implements the Armijo condition for sufficient
% decrease:
%
%           f(x + t*p) < f(x) + c1*t*g(x)'*p
%
% Where t is the step size, p is the search direction, and f(x) and g(x)
% yield the function value and the gradient at x. The implementation of the
% Armijo condition is completed using backtracking. Where t is decrease by
% c2 (where 0<c2<1) until the inequality is satisfied. The default values
% of c1 and c2 are 1E-4 and 0.8. These constitute a 'fine' line-search that
% is happy with small (but sufficient) decreases.
function [t k] = backtrack(fcn,x,p,f,g,setup)
t = 1;% initiate the line search at a pure newton step
if p'*g >= 0
    warning('LineSearch Error: p (p''g = %e)is not a descent direction',p'*g);
end
k = 0;
% Implement The Armijo condition for sufficient decrease
% f(x + t*p) < f(x) + c1*t*g(x)'*p
while fcn(x + t * p) > f + setup.c1*t*(g'*p)
    t = setup.c2*t;% decrease step size by c2 (where 0<c2<1)
    k = k + 1;
%     if k > setup.maxLSits
%         warning('Maxits used in backtracking');
%     end
end
fprintf('\t%1.6f\t%i',t,k);


end



%% Lm-BFGS
%
% Written By Eldad Haber
%
% s(k) = S(:,khat)
% y(k) = Y(:,khat)
%
% where
% n    = #elems,
% len  = min(n , bufferlen),
% khat = mod(n-1-len+k,len)+1,
%
% for k = 1,2,...,len
%
% updating the buffers S,Y works as follows:
%
% updating the buffers with s(n+1) , y(n+1):
% 1) calc top index,     ktop       = mod(n,len)+1
% 2) enter vals,         S(:,ktop)  = s(n+1)
% 3) enter vals,         Y(:,ktop)  = y(n+1)
% 4) increase counter,   n          = n+1
function[d] = bfgs(n,S,Y,H0,d)
len     = min([size(S,2),n]);
ktop    = len;
d       = bfgsrec(ktop,n,len,S,Y,H0,d);
end
function [d] = bfgsrec(k,n,len,S,Y,H0,d)
% BFGS recursion
pcgtol     = 1e-8;
pcgmaxiter = 50;
if k < 1,
    if isnumeric(H0)
        %% solve with pcg
        % M        = spdiags(diag(H0),0,size(H0,1),size(H0,2));
        %H0(1,1) = H0(1,1)+0.2;
        % use sgs precond
%         d = H0\d;
        if length(d) < 70^2
            d = H0\d;
        else
            LL = tril(H0,-1);
            UU = triu(H0,1);
            DD = spdiags(diag(H0),0,size(H0,1),size(H0,2));
            M  = @(f) (UU+DD)\(DD*((LL+DD)\f));
            % solve
            [d] = pcg(H0,d,pcgtol,pcgmaxiter,M);
        end
    else %% assume user supplied solver
        [d,flag] = H0(d);
    end
else
    khat    = mod(n-1-len+k,len)+1;
    gamma   = (S(:,khat)'*d)/(Y(:,khat)'*S(:,khat));
    d       = d - gamma*Y(:,khat);
    d       = bfgsrec(k-1,n,len,S,Y,H0,d);
    d       = d + (gamma - (Y(:,khat)'*d)/(Y(:,khat)'*S(:,khat)))*S(:,khat);
end
end

% A setup structure can be passed in to switch between internal methods.
% The setup structure has fields:
%
% SETUP:
%		c1:              1.000000e-04
%		c2:              0.6000
%		maxStepSize:	 0.1000
%		maxLSits:		 15
%       nbfgs:           15
%       plotIt:          false
%
function setup = ensureSetup(setup)
% This function ensures the correct setup of your structure

names = fieldnames(setup);
trueNames = {'c1','c2','maxStepSize','minStepSize','maxLSits','nbfgs','plotIt'};
for i =1:length(names)
	if ~any(strcmp(names{i},trueNames))
		warning('Setup:FieldNR','Field not recognized: ''%s''',names{i});
	end
end

if ~isfield(setup, 'c1')
	setup.c1 = 1.000000e-04;
end

if ~isfield(setup, 'c2')
	setup.c2 = 0.6000;
end

if ~isfield(setup, 'maxStepSize')
	setup.maxStepSize = 0.1000;
end

if ~isfield(setup, 'minStepSize')
	setup.minStepSize = 1E-4;
end

if ~isfield(setup, 'maxLSits')
	setup.maxLSits = 20;
end

if ~isfield(setup, 'nbfgs')
    setup.nbfgs = 15;
end
if ~isfield(setup, 'plotIt')
    setup.plotIt = false;
end

end
