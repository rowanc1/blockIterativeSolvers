%% Surface Survey design creation in 3D
% Create a simple survey design that is in 3D


% P will be the recierver locations
% Q will be the source locations

% Note on the dimensions of the model m = [x y z] in a left handed
% coordinate system (z goes down!)

function [P Q W] = createSurvey3D(para,design)

design = ensureDesign(para,design);

%% Which survey would you like to do?
% where the full survey takes all permutations of vertical and horizontal
% measurments on the grid for each source pair. Note that the actual
% measurments will be less than this because you cannot measure a voltage
% where you are injecting current. These will be taken care of in a W
% matrix that masks the full data set down into the right dataset.
types = {'perms','monopole'};
switch design.type
    case types{1}%The full permutations of source/reciever locations
        [P Q] = getPerm(para,design);
    case types{2}%monopole
        Q = zeros(para.dims);
        Q(round(para.dims(1)/2),round(para.dims(2)/2),1) = 1;
        Q = Q(:);
        P = Q';%no measurements
    otherwise
        disp('The survey design type was not recognized, choose one of:')
        disp(types);
        error('');
end

W = getW(P,Q);
P = sparse(P);%because n dimensional indexing doesn't work on sparse matricies...
Q = sparse(Q);
Q = Q.*prod(para.dims);%normalize by volume
end

function design = ensureDesign(para,design)
%this function ensures the correct setup of your 'design' structure
if ~isfield(design, 'type')
    disp('Doing the default grid permutation survey');
    design.type = 'perms';
end
if ~isfield(design, 'padding')
    %pad the sides by 15% so work in the middle 70% of your grid
    disp('Adding 15% padding by default');
    design.padding = round(para.dims.*0.15);
end
if ~isfield(design, 'skip')
    %How many cells between electrodes
    disp('Default skipping space used is 2');
    design.skip = 2;
end
if ~isfield(design, 'plotIt')
    design.plotIt = true;
end
end

function [P Q] = getPerm(para,design)
Pbase = zeros(para.dims);
P = [];Q = [];
cnt = 1;
locs = zeros(para.dims(1:2));
% locs(1+design.padding(1):design.skip:end-design.padding(1),[4 6]) = 1;
locs(1+design.padding(1):design.skip:end-design.padding(1),...
     1+design.padding(2):design.skip:end-design.padding(2)) = 1;
if design.plotIt
    figure;
    imagesc(locs);
    title('Surface Electrode Layout')
end
inds = find(locs);
for i = 1:length(inds)-1
    for j = i+1:length(inds)
        Pbase = Pbase.*0;
        [ii ji] = ind2sub(para.dims(1:2),inds(i));
        [ij jj] = ind2sub(para.dims(1:2),inds(j));
        Pbase(ii,ji,1) =  1;
        Pbase(ij,jj,1) = -1;
        P(cnt,:) = Pbase(:)';
        Q(:,cnt) = Pbase(:);
        cnt = cnt +1;
    end
end
end

%% Get W
% W is a masking matrix which, for each RHS (or column of Q), finds the
% locations of the sources that are turned on for the measurment, and finds
% ALL of the recievers that are in conflict with those specific sources
% (i.e. they are in the same model cell).
%
% W is 1 where data is collected and 0 elsewhere.
function [W] = getW(P,Q)
W = ones(size(P,1),size(Q,2));
for i = 1:size(Q,2)%for each RHS
    src = find(Q(:,i));%find which sources are on.
    [srcrec,temp] = find(P(:,src));%if there are any receivers at those locations
    %make sure srcrec only refers to the row index...
    W(srcrec,i) = 0;%turn them off.
end
end


function [P Q W] = getBox2(m,padding)
Pbase = zeros(m);
P = zeros((m(1)-1)*(m(2)-1)*2,prod(m));
cnt = 1;
for i = 1:m(1)-1
    for j = 1:m(2)-1
        Pbase = Pbase.*0;
        Pbase(i, j ,1) =  1;
        Pbase(i,j+1,1) = -1;
        P(cnt,:) = Pbase(:)';
        cnt = cnt +1;
        Pbase = Pbase.*0;
        Pbase( i ,j,1) =  1;
        Pbase(i+1,j,1) = -1;
        P(cnt,:) = Pbase(:)';
        cnt = cnt +1;
        if(j == m(2) - 1);
            Pbase = Pbase.*0;
            Pbase( i ,j+1,1) =  1;
            Pbase(i+1,j+1,1) = -1;
            P(cnt,:) = Pbase(:)';
            cnt = cnt +1;
        end
        if(i == m(1) - 1);
            Pbase = Pbase.*0;
            Pbase(i+1, j ,1) =  1;
            Pbase(i+1,j+1,1) = -1;
            P(cnt,:) = Pbase(:)';
            cnt = cnt +1;
        end
        
    end
end

end
