%% Driver for Solving Multiple RHS
% All work is done on a para file that contains all the information for a
% DC resistivity inversion. Here, para.A is A and para.Q is B and has
% multiple RHSs
load finalPara


%% Use Block CG for n systems
n = 10;
tic;x = bcg(para.A,para.Q(:,1:n),struct('maxits',3000,'tol',1E-5));toc


%% Could compare that to a simple CG algorithm

tic;
for i = 1:n
    x = cg(para.A,para.Q(:,i),struct('maxits',3000,'tol',1E-5));
end
toc


%% Could also do Block Minimum Residual
% With comparable results
tic;x = BlMRes(para.A,para.Q(:,1:n),struct('maxits',3000,'tol',1E-5));toc


%% or Compare that to Minimum Residual
% Uses Ruhe varient on block lanczos, and does not implement block
% deflation inside the algorithm.
% Only thing different here is storage requirements!
tic;
for i = 1:n
    x = MRes(para.A,para.Q(:,i),struct('maxits',3000,'tol',1E-5));
end
toc

%% If you are really bold: Solve it with Block-GMRES!
%With householder reflections
tic;x = blockGMRES(para.A,para.Q(:,1:n),struct('maxits',3000,'tol',1E-5));toc
%Or with Givens rotations!
% Note that there are some problems with this code (works for symmetric
% systems at the moment, not very helpful). But with a few tweaks it should
% work!
tic;x = blockGMRESgiven(para.A,para.Q(:,1:n),struct('maxits',3000,'tol',1E-5));toc


%% The experiments ran looked like this:
% But they take a while to run...

load finalPara
rB = randn(size(para.Q));
for i = 1:size(para.Q,2)
    tic;
    [x Sindiv(i).k Sindiv(i).r Sindiv(i).nrm] = ...
        bcg(para.A,para.Q(:,i),struct('maxits',2000,'tol',1E-5));
    Sindiv(i).time = toc;

    tic;
    [x Sgroup(i).k Sgroup(i).r Sgroup(i).nrm] = ...
        bcg(para.A,para.Q(:,1:i),struct('maxits',2000,'tol',1E-5));
    Sgroup(i).time = toc;
    
    tic;
    [x SindivRand(i).k SindivRand(i).r SindivRand(i).nrm] = ...
        bcg(para.A,rB(:,i),struct('maxits',2000,'tol',1E-5));
    SindivRand(i).time = toc;

    tic;
    [x SgroupRand(i).k SgroupRand(i).r SgroupRand(i).nrm] = ...
        bcg(para.A,rB(:,1:i),struct('maxits',2000,'tol',1E-5));
    SgroupRand(i).time = toc;
    
end



