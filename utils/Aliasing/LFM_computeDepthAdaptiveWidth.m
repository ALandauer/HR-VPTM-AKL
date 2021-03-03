% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function widths = LFM_computeDepthAdaptiveWidth(Camera, Resolution)

%% compute the depth dependent witdth of the anti-aliasing filters
dz = Resolution.depths;
zobj = Camera.fobj - dz;

% avoid division by zero
for i = 1:length(zobj)
    if(zobj(i) == Camera.fobj || zobj(i) == Camera.dof)
        zobj(i) = zobj(i) + 0.00001*Camera.fobj;
    end
end

% z1 -> where the objective will focus
z1 = (zobj * Camera.fobj)./(zobj - Camera.fobj);

% effective radius of the tube lens
tubeRad = Camera.objRad * Camera.Delta_ot*abs(1./z1 - 1/Camera.Delta_ot);

% z2 -> where the tl will focus
z2 = Camera.ftl*(Camera.Delta_ot - z1)./(Camera.Delta_ot - z1 - Camera.ftl);

% main blur (at the mla)
B = tubeRad * Camera.tube2mla .* abs(1./z2 - 1/Camera.tube2mla);

% z3 -> where the mla focuses
z3 = Camera.fm(1)*(Camera.tube2mla - z2)./(Camera.tube2mla - z2 - Camera.fm(1));

% miclolens blur radius
b = Camera.lensPitch/2  * abs(1./z3 - 1/Camera.mla2sensor);
% b = Camera.lensPitch/2  * abs((z3 - Camera.mla2sensor)./z3);

% microlens array to sensor magnification
lambda = z2*Camera.mla2sensor./(Camera.tube2mla * abs(Camera.tube2mla - z2));

% cut-off freq
d = Camera.lensPitch;
f0 = 1/(2*d);

pinhole_filt_rad = d*(abs(lambda)); %1./(2*f0*abs(lambda));
final_rad = abs(pinhole_filt_rad - b);

% Size of filters
widths = min(d/2,final_rad);

%% filter size in object space (voxels)
widthsX = widths * Resolution.TexNnum(2)/d;
widthsY = widths * Resolution.TexNnum(1)/d;

widths = zeros(length(widths), 2);
widths(:,1) = floor(widthsY.*2);
widths(:,2) = floor(widthsX.*2);
widths(mod(widths,2) == 0) = widths(mod(widths,2) == 0) + 1;