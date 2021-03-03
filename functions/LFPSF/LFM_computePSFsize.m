% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function PSFsize = LFM_computePSFsize(maxDepth, Camera)

%% geometric blur radius at the MLA
maxDepth = maxDepth - Camera.offsetFobj;
zobj = Camera.fobj - maxDepth;
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
BlurRad = tubeRad * Camera.tube2mla .* abs(1./z2 - 1/Camera.tube2mla);
PSFsize = ceil(BlurRad/Camera.lensPitch) + 2; %% allow for some extra extent as the psf decays smoothly compared to rays prediction
disp(['Size of PSF radius ~= ' num2str(PSFsize) ' [microlens pitch]' ]);