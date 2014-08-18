%% Create Data (3D)
% Adds random noise to the data.
% The amount of noise added is taken from design.error

function [dataOut] = createData3D(data,design)


outlierError = 0.25;
if design.outliers > 0
    f = find(data);
    p = randperm(nnz(data));
    for i = 1:design.outliers
        data(f(p(i))) = data(f(p(i))) + randn(1).*(outlierError.*data(f(p(i))));
    end
end
design = ensureDesign(design);
types = {'normal'};
switch design.type
    case types{1}%normal
        noise = randn(size(data)).*(design.error.*data);
        dataOut = data + noise;
    otherwise
        disp('The data design type was not recognized, choose one of:')
        disp(types);
        error('');
end

end

function design = ensureDesign(design)
%this function ensures the correct setup of your 'design' structure
if ~isfield(design, 'error')
    disp('Adding 2% noise to the data.');
    design.error = 0.02;
end
if ~isfield(design, 'type')
    disp('Adding random *normal* noise by default.');
    design.type = 'normal';
end
if ~isfield(design, 'outliers')
    % Not adding any outling data...
    design.outliers = 0;
end
end

