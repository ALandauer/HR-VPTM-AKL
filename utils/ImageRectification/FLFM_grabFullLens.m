% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu 

function MLAmask = FLFM_grabFullLens(Camera, Resolution)

% this builds a mask for the MLA which only keeps these elemental images that
% fully fit into the sensor

circ = @(x,y,r) (x.^2.+y.^2.<r.^2)*1.0;
maskSensor = zeros(Resolution.Nnum);
y = linspace(-floor(size(maskSensor,1)/2), floor(size(maskSensor,1)/2), size(maskSensor,1));
x = linspace(-floor(size(maskSensor,2)/2), floor(size(maskSensor,2)/2), size(maskSensor,2)); 
[X,Y] = meshgrid(x,y);
maskSensor = circ(X, Y, Camera.spacingPixels/2);

ylength = Resolution.sensorSize(1);
xlength = Resolution.sensorSize(2);

% centers
LensletCenters(:,:,1) = round(Resolution.LensletCenters(:,:,2));
LensletCenters(:,:,2) = round(Resolution.LensletCenters(:,:,1));

% activate lenslet centers outside the sensor size
MLcenters = zeros(ylength, xlength);
for a = 1:size(LensletCenters,1) 
    for b = 1:size(LensletCenters,2) 
        if( (LensletCenters(a,b,1) < Resolution.Nnum(1)/2) || (LensletCenters(a,b,1) > ylength - Resolution.Nnum(1)/2) || ...
              (LensletCenters(a,b,2)< Resolution.Nnum(2)/2) || (LensletCenters(a,b,2) > xlength - Resolution.Nnum(2)/2) )
         MLcenters( LensletCenters(a,b,1), LensletCenters(a,b,2)) = 1;
        end
    end
end
% apply the mask to the centers outside the sensor size
MLAmask  = conv2(MLcenters(1:ylength, 1:xlength), maskSensor, 'same');

% invert the mask
MLAmask = ~MLAmask;