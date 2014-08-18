%% ParaClass
% This structure holds all of the parameters about your problem and
% the majority of the code to run that problem.
%
% Initiate the class by calling:
%    para = ParaClass([dim1 dim2 dim3]);
% Where dim# is the dimension of each axis.
%
% Next you can create a survey design
%    para = para.createSurvey(struct());
% Notice that you MUST rewrite para when you are adding fields to it.
%
% See Also createData3D createSurvey3D createModel3D
classdef ParaClass
    properties
        dims%Holds the dimensions of the model. Assumes regular grid from 0-1 in each direction
        D%Holds the observed data
        P%Reciver locations
        Q%Source locations normalized by volume
        W%Mask for where you actually collect data.
        DIV%Divergence operator;
        GRAD%Gradient operator
        GRADw% weighted gradient with smallness
        Av%Cell centered averaging matrix.
        
        % Hold some information about the current model
        m% model
        mref% reference model
        A%Holds the LHS for the current model
        % Regularization parameters
        beta
        
        objFnc %chooses the objective function that you use.
    end
    
    properties (SetAccess = private)
        %% Quick Functions
        % These are defined when you first call para.
        ddx%the dirivative in one dimension
        ddxc%cell centered gradient in one dimension, includes nuemann boundary conditions
        av%The averaging matrix
        kron3%do a quick 3D kron-product
        sdiag%short form for 1 sparse diagonal
        
    end
    methods
        function para = ParaClass(dims_in)
            
            if nargin > 0 && length(dims_in) ==3
                para.dims = dims_in;
            else
                error('You must include a 3 element dimension vector')
            end
            
            % Create all of the quick function calls
            para.ddx   = @(n)(n*spdiags(ones(n+1,1)*[-1,1],[0,1],n,n+1));%this assumes 0:1 grid in all dims
            % ddxc puts in nuemann boundary conditions.
            % this is actually the negative of the gradient
            para.ddxc  = @(n)(n*spdiags([[-ones(n-1,1);0],[0;ones(n-1,1)]],[-1,0],n+1,n));
            para.av    = @(n)(spdiags(ones(n+1,1)*[1,1]/2,[-1,0],n+1,n));
            para.kron3 = @(A,B,C) kron(A,kron(B,C));
            para.sdiag = @(a) spdiags(a(:),0,numel(a),numel(a));
            
            
            % Create DIV and Av
            dm = para.dims;%short form for dimensions
            para.DIV = [para.kron3(speye(dm(3)),speye(dm(2)),para.ddx(dm(1))), ...
                para.kron3(speye(dm(3)),para.ddx(dm(2)),speye(dm(1))), ...
                para.kron3(para.ddx(dm(3)),speye(dm(2)),speye(dm(1)))];
            
            para.GRAD = [para.kron3(speye(dm(3)),speye(dm(2)),para.ddxc(dm(1))); ...
                para.kron3(speye(dm(3)),para.ddxc(dm(2)),speye(dm(1))); ...
                para.kron3(para.ddxc(dm(3)),speye(dm(2)),speye(dm(1)))];

            para.GRADw = para.createGRADw(1,1,0.1,0.1,0);%less in the z direction.
            
            % Averaging matrix from cell centers to nodes.
            para.Av  = [para.kron3(speye(dm(3)),speye(dm(2)),para.av(dm(1))); ...
                para.kron3(speye(dm(3)),para.av(dm(2)),speye(dm(1))); ...
                para.kron3(para.av(dm(3)),speye(dm(2)),speye(dm(1)))];
            
            para.objFnc = 'l2';
        end
        % Weighted Gradient in XYZ and Smallness wrt mref
        % where d is the depth weighting exponent in 1/(z_0+z_ij)^d
        % Set d=0 to have a straight smallness term
        function Gw = createGRADw(para,x,y,z,s,d)
            dm = para.dims;%short form for dimensions
            z0 = 1;
            depth = repmat(reshape(linspace(0,1,dm(3)),[1 1 dm(3)]),[dm(1) dm(2) 1]);
            Gw = [x*para.kron3(speye(dm(3)),speye(dm(2)),para.ddxc(dm(1))); ...
                y*para.kron3(speye(dm(3)),para.ddxc(dm(2)),speye(dm(1))); ...
                z*para.kron3(para.ddxc(dm(3)),speye(dm(2)),speye(dm(1))); ...
                s*para.sdiag((z0+depth(:)).^(-d))];%smallness
        end
        %% Create Survey
        % Creates a new survey design with the design parameters in a
        % structure field.
        % See createSurvey3D for more information on the structure field.
        function para = createSurvey(para,surveyDesign)
            if nargin == 2 && isstruct(surveyDesign);
                [para.P,para.Q,para.W] = createSurvey3D(para,surveyDesign);
            else
                help('createSurvey3D')
                error('You need to input the surveyDesign.')
            end
        end
        %% Add Survey Info
        % Integrate P, Q, and W into para.
        function para = addSurveyInfo(para,p,q,w)
            if nargin == 4
                para.P = p;
                para.Q = q;
                para.W = w;
            else
                error('You need to input the P Q and W mats');
            end
        end
        %% Create Model
        % This will make a new model based on the fields of the design
        % structure.
        % See Also createModel3D
        function [para] = createModel(para,design)
            if nargin == 2 && isstruct(design)
                new_m = createModel3D(para,design);
                para = newModel(para,new_m);
            else
                help('createModel3D')
                error('You must input a design structure to make a model. If you have one in mind already, you can call para.newModel(m)');
            end
        end
        %% Create Data
        % This will create data and add noise based on the design
        % parameters.
        % See createData3D for more information.
        function [para] = createData(para,design,sd,e)
            if nargin > 1 && isstruct(design)
                para.D = para.W.*0;%zeros out the data so that misfit is true data.
                para.W(para.W ~= 0) = 1;%Reset the weighting matrix to be 0s and 1s
                b = para.beta;
                para.beta = 0;%set regularization to zero when you are creating data.
                [temp,temp,data] = para.getMisfit();
                para.beta = b;%redefine the regularization parameters
                para.D = createData3D(data,design);
                if nargin > 2
                    para.W(para.W ~= 0) = 1./(abs(para.D(para.W ~= 0))*sd+e);
                else
                    warning('No data weighting was applied')
                end
            else
                help('createData3D');
                error('You must input a design field, that will go into createData3D.');
            end
            
        end
        %% New Model
        % Given a new model (m) this code will calculate the new A(m)
        function [para] = newModel(para,m)
            if nargin == 2 && numel(m) == prod(para.dims)
                para.m = m(:);
                para.A = para.getA(m);
            else
                error('You must input a model (m) of the right size');
            end
        end
        function [a] = getA(para,m,fixNull)
            %% TODO??
            % Make the inversion invariant to working in resistivity or
            % conductivity by working in log space. This will incur a
            % few new terms in the grad equation because of chain rule!
            % sig = e^m;
            a = -para.DIV* para.sdiag(para.Av*exp(m(:))) * para.GRAD;
            
            % Fix the nullspace of the matrix by adding 1 in the first element
            if nargin < 3 || fixNull% pass in para.getA(m,false) to turn off null fix.
                a(1) = a(1) + (sign(a(1))*1);%make sure it is the same sign.
            end
        end
        %% Get Misfit
        % Returns the data misfit for the current model
        function [misfit,dmisfit,R,Y] = getMisfit(para,m)
            if nargin > 1%if you give me a new model here, update para, and recalculate A
                para = para.newModel(m);
            end
            
            
            
            
            %Y is the potential field everywhere in the model space
            Y      = para.A\para.Q;%Solve for all RHS at once.
            R      = para.P*Y - para.D;%Subtract off the data.
            R      = para.W.*R;%Zero out the ones we don't care about.
            
            switch para.objFnc
                case 'l2'%least squares
                    misfit = 0.5*(R(:)'*R(:));
                case 'l1'%L1 norm approximation
                    eps = 1E-5;
                    misfit = sum(sqrt(R(:).^2+eps));
                case 'st'% Student T distribution
                    w = 1;
                    misfit = sum(log(w+R(:).^2));
                otherwise
                    error('Objective function not supported. Choose: ''l2'', ''l1'' or ''st''.');
            end
            if nargout > 1% compute derivatives if needed
                % get the objective gradient
                dmisfit = 0;
                for i=1:size(para.Q,2)%Solve for dirivative of each RHS
                    %You just include W' in the dirivative. where W can be
                    %thought of as a diagonal matrix, and it will flip to
                    %the other side and you can just multiply it here!
                    switch para.objFnc
                        case 'l2'%least squares
                            dmisfit = dmisfit + para.Jt(Y(:,i),...
                                para.W(:,i).*R(:,i));
                        case 'l1'%L1 norm approximation
                            dmisfit = dmisfit + para.Jt(Y(:,i),...
                                para.W(:,i).*(1./sqrt(R(:,i).^2+eps)).*R(:,i));
                        case 'st'% Student T distribution
                            dmisfit = dmisfit + para.Jt(Y(:,i),...
                                para.W(:,i).*(2./(w+R(:,i).^2)).*R(:,i));
                        otherwise
                            error('Objective function not supported. Choose: ''l2'', ''l1'' or ''st''.');
                    end
                end
            end
%             fprintf('The misfit is: %e \n',misfit);
            
            
            % Add regularization
            [regF regG] = para.getRegMisfit(para.m);
            misfit = misfit + regF;
            if nargout > 1
                dmisfit = dmisfit + regG;
            end
            
        end
        function [f,g] = getRegMisfit(p,m)
            if isempty(p.beta)
                warning('Regularization parameter not set. para.beta = 1E-3;');
                p.beta = 1E-3;
            end
            if isempty(p.mref)
                warning('Reference Model is not set. para.mref = 0;');
                p.mref = 0;
            end
            f = p.beta/2*(norm(p.GRADw*(m-p.mref))^2);
            g = p.beta.*p.GRADw'*(p.GRADw*(m-p.mref));
        end
        %% Jacobian Functions
        % This function updates the sensitivity matricies (of which there
        % are three).
        % Where u == A\Q(:,i)
        % and v is the multiplying vector.
        % To get the full gradient, use para.dCdm(u,v,m) where m is the
        % model.
        function [Jv] = J(p,u,v)
            Jv = -p.P*(p.A\p.dCdm(u,v,p.m));
        end
        function [Gv] = dCdm(p,u,v,m)
            Gv = -(p.DIV*(p.sdiag(p.GRAD*u) * (p.Av*(p.sdiag(exp(m))*v))));
        end
        function [Jtw] = Jt(p,u,w)
            y = ((p.A')\(p.P'*w));
            Jtw = -p.sdiag(exp(p.m))*((p.Av')*(p.sdiag(p.GRAD*u) * ((p.GRAD)*y)));
        end
        
    end
end % classdef
