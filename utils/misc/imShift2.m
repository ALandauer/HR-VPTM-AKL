% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2019 Anca Stefanoiu and Josue Page

function newImg = imShift2(Img, ShiftX, SHiftY)

eqtol = 1e-10;

xlength = size(Img,1);
ylength = size(Img,2);

if abs(mod(ShiftX,1)) > eqtol || abs(mod(SHiftY,1)) > eqtol
   error('SHIFTX and SHIFTY must be integers');
end
    
ShiftX = round(ShiftX);
SHiftY = round(SHiftY);

newImg = zeros(xlength, ylength, size(Img,3) );

if ShiftX >=0 && SHiftY >= 0
    newImg( (1+ShiftX:end), (1+SHiftY:end),:) = Img( (1:end-ShiftX), (1:end-SHiftY),:);
elseif ShiftX >=0 && SHiftY < 0
    newImg( (1+ShiftX:end), (1:end+SHiftY),:) = Img( (1:end-ShiftX), (-SHiftY+1:end),:);
elseif ShiftX <0 && SHiftY >= 0
    newImg( (1:end+ShiftX), (1+SHiftY:end),:) = Img( (-ShiftX+1:end), (1:end-SHiftY),:);
else
    newImg( (1:end+ShiftX), (1:end+SHiftY),:) = Img( (-ShiftX+1:end), (-SHiftY+1:end),:);
end