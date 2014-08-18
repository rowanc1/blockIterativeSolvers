%% Create a model
% Has various parameters that you can adjust in the design call.

function [m] = createModel3D(para,design)

design = ensureDesign(design);
[X Y Z] = ndgrid(linspace(0,1,para.dims(1)),linspace(0,1,para.dims(2)),linspace(0,1,para.dims(3)));
m = zeros(para.dims);%initialize the model

%types = {'homo','block','horzLayers','vertLayers'};
types = {'homo','block','eh','vertLayers'};

switch design.type
    case types{1}%Homogeneous
        m(:) = design.values(1);
    case types{2}%block
        m = createBlock(para,design);
    case types{3} % 2 blocls
        nn = para.dims;
        m = ones(nn)*log(1e-2);
        m(fix(nn(1)/4):fix(nn(1)/2),fix(nn(2)/4):fix(nn(2)/2),3:5) = 0;
        m(fix(nn(1)/2):fix(3*nn(1)/4),fix(nn(2)/2):fix(3*nn(2)/4),2:4) = log(1e-4);
        m = m(:);
    otherwise
        disp('The model type was not recognized, choose one of:')
        disp(types);
        error('');
end
if design.addError
    disp('Adding error to the model.');
    m = createData3D(m,design.error);
end
if design.plotIt
    montageArray(reshape(m,para.dims));
end
m = m(:);%Unwrap it.
end

function design = ensureDesign(design)
%this function ensures the correct setup of your 'design' structure


if ~isfield(design, 'type')
    disp('Creating a homogeneous model');
    design.type = 'homo';
end
if ~isfield(design, 'size')
    design.size = [3,3,3];
end
if ~isfield(design, 'position')
    design.position = 'centered';%or random
end
if ~isfield(design, 'values')
    design.values = 1E-5;
end
if ~isfield(design, 'addError')
    design.addError = false;
end
if ~isfield(design, 'error')
    design.error = struct();
end

if ~isfield(design, 'plotIt')
    design.plotIt = true;
end

end


function m = createBlock(para,design)
m = ones(para.dims);
switch design.position
    case 'center'
        d = round((para.dims-(design.size-1))./2);
        d(3) = 1;
    case 'middle'
        d = round((para.dims-(design.size-1))./2);
    case 'corner'
        d = [1 1 1];
    case 'random'
        d = floor((para.dims-(design.size-1)).*rand(1,3));
        d(d==0) = 1;
    otherwise
        error('The position of the block can either be centered or random');
end
m(d(1):(d(1)+design.size(1)),...
  d(2):(d(2)+design.size(2)),...
  d(3):(d(3)+design.size(3))) = 0;
m(m==1) = design.values(1);
m(m==0) = design.values(2);

end


