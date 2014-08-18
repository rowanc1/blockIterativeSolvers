function testCode(drvIt)
types = {'forward','gradient','adjoint','all'};
if ~exist('drvIt','var')
    disp(types)
    type = types{input('Choose a test type (1-4):  ')};
else
    type = types{drvIt(1)};
end
switch type
    case types{1}% test the forward
        disp('*****     FORWARD MODEL TEST     *****');
        fprintf('    h   \t ||err||inf \n')
        for i = [8, 16, 32, 64]%, 128]
            
            h = 1/i;
            t = h/2:h:1-h/2;
            [x,y,z] = ndgrid(t,t,t);
            
            para = ParaClass([i i i]);
            
            A = para.getA(zeros(i^3,1),false);%zeros because exp(m), don't apply a nullSpace fix to A
            
            u = cos(2*pi*x).*cos(2*pi*y).*cos(2*pi*z);%something with nuemann BC
            
            b = 3*(2*pi)^2*u;%Note non negative because grad is really -grad fix this.
            
            bn = A*u(:);%calculate numerical RHS
            
            err = bn - b(:);%Compare to true RHS
%             montageArray(reshape(err,i,i,i)); colorbar
%             pause(0.5);figure(gcf)
            fprintf('    %3.2e \t %3.2e\n',h,norm(err,'inf'));
        end
        
    case types{2}
        %% Gradient Test
        % The Gradient test, tests whether adding gradient information to
        % the problem decreases the error quadratically.
        %
        % Where
        % |c(u+hv,m)-c(u, v)|
        % |c(u+hv,m)?c(m,u)?h*G_u*c*v|
        
        para = ParaClass([8 8 8]);
        para.mref = 0;
        objFcnTypes = {'l2','l1','st'};
        if exist('drvIt','var') && length(drvIt)==2
            objType = objFcnTypes{drvIt(2)};
        else
            disp(objFcnTypes);
            objType = objFcnTypes{input('Choose and Objective Function (1-3):  ')};
        end
        para.objFnc = objType;
        %Create some random survey.
        surveyDesign = struct('type','perms','padding',[1 1],'skip',2,'plotIt',0);
        para = para.createSurvey(surveyDesign);
        [n1,n2] = size(para.W);
        para.W = rand(n1,n2);
        para.beta = 0;%look at this without regularization
        
        m    = randn(8^3,1);
        para = para.newModel(m);
        para.D = 0;%just zero for data
        
        [f,df]  = para.getMisfit(m);
        dm   = randn(8^3,1);
        
        fprintf('*********     GRADIENT TEST - %s     *********\n',objType);
        fprintf('    h   \t   fp-f   \tfp-f-h*G*v  \n')
        for i = 1:10
            h = 10^(-i);
            f1 = para.getMisfit(m+h*dm);
            fprintf('%3.2e \t %3.2e \t %3.2e\n',h,abs(f1-f),abs(f1-f-h*df'*dm));
        end
        
    case types{3}
        %% Adjoint Test
        % Look for the wtJv and the transpose (vtJtw) such that the
        % difference is at machine error.
        para = ParaClass([8 8 8]);
        
        %Create some random survey.
        surveyDesign = struct('type','perms','padding',[1 1],'skip',2,'plotIt',0);
        para = para.createSurvey(surveyDesign);
        para.beta = 0;%look at this without regularization
        m    = randn(8^3,1);
        para = para.newModel(m);
        
        
        %Create some random vectors to test with
        w = randn(size(para.P,1),1);
        v = randn(size(para.Av,2),1);
        u = randn(size(para.A,1),1);
        wtJv  = w'*para.J(u,v);
        vtJtw = v'*para.Jt(u,w);
        disp('*********     ADJOINT TEST     *********');
        fprintf('wtJv - vtJtw = %e\n',abs(wtJv - vtJtw));
    case types{4}
        % Test Everything
        testCode(1)
        testCode([2 1])
        testCode([2 2])
        testCode([2 3])
        testCode(3)
    otherwise
        disp('The test type was not recognized, choose one of:')
        disp(types);
        error('');
end

end
