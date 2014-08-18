%% GMRES(m) - Restarted Generalized Minimum Residual
% The algorithm iterativly solves the system Ax=b to a specified tolerance
% (tol) within the maximum number of iterations specified (maxits). The
% algorithm will restart the algorithm every *restart* times, this can be
% handy if memory is going to be an issue.
%
% The function uses a setup structure with a number of flags, explained
% below. To get this structure, run the function with no inputs. To use
% preconditioning and to create an LU factorization inside this function,
% leave L and U blank, and let precond=true and iluType = p in ilup(A,p).
%
%   setup = 
%
%           maxits: 10          Maximum Iterations
%          restart: 10          Restart Number
%              tol: 1.0000e-06  Solve Tolerance
%            isSym: 0           Flag for A symmetry
%     showComments: 0           Flag for showing comments
%     useBackSlash: 1           Flag for using \ in the preconditioner slve
%          precond: 0           Flag for using preconditioning
%          iluType: 0           p in ILU(p) is preconditioning is used 
%                L: MATRIX      L matrix (can be left blank)
%                U: MATRIX      U matrix (can be left blank)
%
% If the system is symmetric use the boolean isSym flag to change from
% solving arnoldi to solving lanczos.
%
% The boolean flag comments can be used to output the residual at every
% step.
%
% The following code uses Given's rotations to find the residual, and the
% real residual is only calculated outside the inner loop.
%
%
% Based on class notes - CPSC 517, and Saad - Iterative Methods

function [X] = blockGMRESgiven(A,B,setup)
%% Set up the arguments

if nargin <2
    disp('You must include at least 2 arguments: A and b');
    X = ensureSetup(struct());
    return;
elseif nargin <3
    setup = ensureSetup(struct());
else
    setup = ensureSetup(setup,A);
end

% How many RHSs?
p = size(B,2);
%% GMRES(restart)
% Set up Arguments
Xk = B.*0;%start with 0 as initial guess
s = zeros(setup.restart,p);c = s;%Givens rotations
R = zeros(setup.restart+1,setup.restart);%upper triangular matrix in norm(g-R*yk)
G = zeros(setup.restart+1,p);
Q = zeros(size(B,1),setup.restart+1);%set up Q and H for Arnoldi
H = zeros(setup.restart+1,setup.restart);
log_nrm_rk = inf;
logTol = log10(setup.tol);
cnt = 1;
loops = ceil(setup.maxits/setup.restart);
%% The outer loop
for loop = 1:loops
    %setup initial guess
    X0 = Xk;
    % Setup Arnoldi
    if setup.precond%left preconditioning!
        %if you are using preconditioning, do a quick forward and back
        %solve to find your R_0* = M^-1 * R_0
        [q r] = qr(setup.U\(setup.L\(B-A*X0)),0);
    else
        [q r] = qr(B-A*X0,0);
    end
    Q(:,1:p) = q;
    G = G.*0;
    G(1:p,1:p) =  r;
    if setup.showComments;fprintf('Restarting. norm of resid: %e\n',G(1));end
    %% The Inner Loop
    for k = p:setup.restart
        j = k-p+1;
        if setup.precond
            z = setup.U\(setup.L\(A*Q(:,j)));
        else
            z = A*Q(:,j);
        end
        vals = max(1,k-2*p):k;
        for i = vals
            H(i,j) = Q(:,i)'*z;
            z = z - H(i,j).*Q(:,i);
        end
        H(k+1,j) = norm(z);
        if abs(H(k+1,j)) < setup.tol;
            fprintf('Terminated at k = %i, H(%i,%i) == 0\n',j,k+1,j);
            % Don't break until you have calculated g and R, the other
            % statments will break out of the loop (just don't divide by
            % zero!)
        end
        Q(:,k+1) = z./H(k+1,j);%find q_(k+1)
        
        % QR Factorize H using Griven's Rotations
        R(1:k+1,j) = H(1:k+1,j);%copy the last column
        % Apply all the previous rotations to the k-th column of H.
        % this is a simplified matrix-vector product Qj'*R_k
        vals = max(1,j-2*p):j-1;
        for rot1 = vals
            for rot2 = p:-1:1
                pls = rot1+rot2;
                mns = rot1+rot2-1;
                Rrj      =  c(rot1,rot2).*R(mns,j) + s(rot1,rot2).*R(pls,j);
                R(pls,j) = -s(rot1,rot2).*R(mns,j) + c(rot1,rot2).*R(pls,j);
                R(mns,j) = Rrj;
            end
        end
        % Solve for the current Givens Rotations
        for rot2 = p:-1:1
            pls = j+rot2;
            mns = j+rot2-1;
            x = R(mns:pls,j);
            r = norm(x);
            c(j,rot2) = x(1)/r;
            s(j,rot2) = x(2)/r;
            % Apply the rotation to the sub-diagonal element
            R(mns:pls,j) = [r;0];
            % Apply the rotation to g
            for col = rot2:p%Due to triangular starting!
                Gmns       =  c(j,rot2).*G(mns,col) + s(j,rot2).*G(pls,col);
                G(pls,col) = -s(j,rot2).*G(mns,col) + c(j,rot2).*G(pls,col);
                G(mns,col) = Gmns;
            end
        end
        % Some testing stuff
%         Yk1 = triu(R(1:j,1:j))\G(1:j,:);
%         Yk2 = H(1:k+1,1:j)\G0(1:k+1,:);
%         n(1) = norm(Yk1(:)-Yk2(:));
%         Xk = X0 + Q(:,1:j)*Yk2(1:j,:);
%         n(2) = norm(B-A*Xk) - norm(G(j+1:end,:));
%         if sum(n) > 1E-3;%setup.tol
%             keyboard
%         end
        
        % Approximate the Norm
        % Only do this if you are close or evey p iterations
        if log_nrm_rk-logTol< 2 || rem(cnt-1,p) == 0
            nrm_rk = norm(G(j+1:end,:));
            log_nrm_rk = log10(nrm_rk);
        end
        if setup.showComments;fprintf('norm(b-A*x) = %e\n',nrm_rk);end
        
        %% Stopping Criteria
        if(nrm_rk < setup.tol)
            break;
        end
        cnt = cnt + 1;
        if cnt > setup.maxits
            break;
        end
    end %inner loop
    %% Solve System for yk (do this inline so you don't redefine variables)
    % Use back substitution to find yk taking triu ensures there aren't any
    % round off errors below the diagonal
    Yk = triu(R(1:j,1:j))\G(1:j,:);
    
    %% Evaluate xk
    Xk = X0 + Q(:,1:j)*Yk;
    if(nrm_rk < setup.tol)
        break;
    end
end %outer loop

%% Output a few things
if(nrm_rk < setup.tol)
    if loop == 1
        fprintf('gmres CONVERGED at iteration %i. \nnorm(b-A*x) = %e\n',k,norm(B-A*Xk));
    else
        fprintf('gmres(%i) CONVERGED at iteration %i(%i). \nnorm(b-A*x) = %e\n',setup.restart,loop,k,norm(B-A*Xk));
    end
else
    if loop == 1
        fprintf('\n\n*****************************\n\ngmres DID NOT converge by iteration %i. \nnorm(b-A*x) = %e\n\n*****************************\n\n',k,norm(B-A*Xk));
    else
        fprintf('\n\n*****************************\n\ngmres(%i) DID NOT converge by iteration %i(%i). \nnorm(b-A*x) = %e\n\n*****************************\n\n',setup.restart,loop,k,norm(B-A*Xk));
    end
end
X = Xk;

end

%% Some Ugly setup stuff Hidden down here...
function setup = ensureSetup(setup,A)
%this function ensures the correct setup of your 'setup' structure
if ~isfield(setup, 'maxits')
    setup.maxits = 10;
end
if ~isfield(setup, 'restart')
    setup.restart = setup.maxits;
end
if ~isfield(setup, 'tol')
    setup.tol = 1E-6;
end
if ~isfield(setup, 'showComments')
    setup.showComments = false;
end
if ~isfield(setup, 'precond')
    if ~isfield(setup, 'L')
        setup.precond = false;
    elseif isfield(setup,'U')% L and U exist
        setup.precond = true;
    else
        setup.precond = false;
    end
    setup.iluType = 0;% just put this in to fill...
elseif setup.precond && ~isfield(setup, 'L') && isfield(setup, 'iluType')
    [setup.L,setup.U] = ilu(A,setup.iluType);%create your L U factorization
elseif setup.precond && isfield(setup, 'L') && isfield(setup, 'U')
    %good to go! - I hope!
    %Could put some checks in here...
else
    setup.precond = false;
end
end
