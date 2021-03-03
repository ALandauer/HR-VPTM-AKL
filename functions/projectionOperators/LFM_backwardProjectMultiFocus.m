% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function BackProjection = LFM_backwardProjectMultiFocus( Ht, projection, lensCenters, Resolution, texSize, range)

lensTypesCount = size(Ht,2); %3
nDepths = size(Ht{1,1},3);

%% perform lensTypesCount backprojections and accumulate the intermediate result
BackProjection = zeros(texSize(1), texSize(2), nDepths);
for i = 1:lensTypesCount
    fprintf(['\nBackProject lens type: ', num2str(i), '/', num2str(lensTypesCount)]);
    lensCurrentType = LFM_extractLensType(lensCenters, i);
    BackProjection = BackProjection + LFM_backwardProject(Ht{1,i}, projection, lensCurrentType, Resolution, texSize, range);
end