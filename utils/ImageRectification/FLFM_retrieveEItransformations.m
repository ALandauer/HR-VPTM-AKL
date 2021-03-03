% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu 

function transformationsStack = FLFM_retrieveEItransformations(LensletGridModel, calibrationImage)

LF = FLFM_extractEI(LensletGridModel, calibrationImage);

% register EIs to central one to fnd the shift between them
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;

fprintf('\nFinding lenslets centers')
transformationsStack = cell(size(LF,3), size(LF,4));
fixed = LF(:,:,ceil(size(LF,3)/2), ceil(size(LF,4)/2));
for i = 1 : size(LF,3)
    for j = 1 : size(LF,4)
        fprintf('.')
        warning('off')
        tform = imregtform(LF(:,:,i,j), fixed, 'translation', optimizer, metric);
        warning('on')
%         LF(:,:,i,j) = imwarp(LF(:,:,i,j),tform,'OutputView',imref2d(size(fixed)));
%         imshowpair(fixed, LF(:,:,i,j),'Scaling','joint')
        transformationsStack{i,j} = tform; % store the transformation to retrieve the EI centers offsets 
    end
end
