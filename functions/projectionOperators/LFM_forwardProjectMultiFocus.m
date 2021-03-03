% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function Projection = LFM_forwardProjectMultiFocus( H, realSpace, lensCenters, Resolution, imgSize, range)

lensTypesCount = size(H,2); %3

%% perform lensTypesCount projections and accumulate the intermediate result
Projection = zeros(imgSize(1), imgSize(2));
for i = 1:lensTypesCount
    fprintf(['\nForwardProject lens type: ', num2str(i), '/', num2str(lensTypesCount)]);
    lensCurrentType = LFM_extractLensType(lensCenters, i);
    Projection = Projection + LFM_forwardProject(H{1,i}, realSpace, lensCurrentType, Resolution, imgSize, range);
end