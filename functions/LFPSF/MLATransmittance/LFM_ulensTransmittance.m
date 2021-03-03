% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function ulensPattern = LFM_ulensTransmittance(Camera, Resolution)

% Compute lens transmittance function for one micro-lens (consitent with Advanced Optics Theory book and Broxton's paper)
ulensPattern = zeros( length(Resolution.yMLspace), length(Resolution.xMLspace), length(Camera.fm) );
for j = 1: length(Camera.fm)
    for a=1:length(Resolution.yMLspace)
        for b=1:length(Resolution.xMLspace)
            x1 = Resolution.yMLspace(a);
            x2 = Resolution.xMLspace(b);
            xL2norm = x1^2 + x2^2;
            ulensPattern(a,b,j) = exp(-1i*Camera.k/(2*Camera.fm(j))*xL2norm);
        end
    end
end

% Mask the pattern, to avoid overlaping when applying it to the whole image 
for j = 1:length(Camera.fm)
    patternML_single = ulensPattern(:,:,j);
    if (Resolution.maskFlag == 1)
        patternML_single(Resolution.sensMask == 0) = 0;
    else
        % circ lens shapes with with blocked light between them
        [x,y] = meshgrid(Resolution.yMLspace, Resolution.xMLspace);
        patternML_single((sqrt(x.*x+y.*y) >= Camera.lensPitch/2  - 3)) = 0; % TODO: hardcoded 3 
    end
    ulensPattern(:,:,j) = patternML_single;
end