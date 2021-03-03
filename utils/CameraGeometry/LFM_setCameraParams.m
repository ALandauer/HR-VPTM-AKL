% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function Camera = LFM_setCameraParams(configFile, newSpacingPx)

%% set LFM parameters:
%%%% objective params
% M-> objective magnification
% NA-> objective aperture
% ftl-> focal length of tube lens (only for microscopes)

%%%% sensor
% lensPitch-> lenslet pitch
% pixelPitch-> sensor pixel pitch

%%%% MLA params
% gridType-> microlens grid type: "reg" -> regular grid array; "hex" -> hexagonal grid array 
% focus-> microlens focus: "single" -> all ulens in the array have the same focal length; 'multi' -> 3 mixed focal lenghts in a hex grid
% fml-> focal length of the lenslets
% uLensMask-> 1 when there is no space between ulenses (rect shape); 0 when there is space between ulenses (circ aperture)

%%%% light characteristics
% n-> refraction index (1 for air)
% wavelenght-> wavelenght of the the emission light

%%%%% distances
% plenoptic-> plenoptic type: "1" for the original LFM configuration (tube2mla = ftl); "2" for defocused LFM (tube2mla != ftl) design
% tube2mla-> distance between tube lens and MLA (= ftl in original LFM design)
% mla2sensor-> distance between MLA and sensor

Camera = ReadYaml(configFile);
if(iscell(Camera.fm))
    Camera.fm = cell2mat(Camera.fm);
end
%% compute extra LFM configuration specific parameters
% for regular grid we only need to compute 1/4 of the psf pattens, due to symmetry
if(strcmp(Camera.gridType, 'reg'))
    Camera.range = 'quarter';
elseif(strcmp(Camera.gridType, 'hex'))
    Camera.range = 'full';
end

Camera.fobj = Camera.ftl/Camera.M;  %% focal length of objective lens
Camera.Delta_ot = Camera.ftl + Camera.fobj; %% obj2tube distance

% ulens spacing = ulens pitch
spacingPx = Camera.lensPitch/Camera.pixelPitch;
if(strcmp(newSpacingPx, 'default'))
    newSpacingPx = spacingPx;
end
Camera.spacingPx = spacingPx;
Camera.newSpacingPx = newSpacingPx; % artificial spacing, usually lower than spacingPx for speed up
Camera.newPixelPitch = (Camera.lensPitch/newSpacingPx);

Camera.k = 2*pi*Camera.n/Camera.WaveLength; %% wave number

objRad = Camera.fobj * Camera.NA; % objective radius

if (Camera.plenoptic == 2)
    
    obj2tube = Camera.fobj + Camera.ftl; %% objective to tl distance
    
    % check if tube2mla is provided by user, otherwise compute it s.t. the F number matching condition is satisfied
    if Camera.tube2mla == 0
        if Camera.mla2sensor == 0
            error('At least one of the "tube2mla" or "mla2sensor" distances need to be provided in the - plenoptic == 2 - case');
        else
            tube2mla = computeTube2MLA(Camera.lensPitch, Camera.mla2sensor, obj2tube, objRad, Camera.ftl);
            Camera.tube2mla = tube2mla;
        end
    end
    
    dot = Camera.ftl * Camera.tube2mla/(Camera.tube2mla-Camera.ftl); %% depth focused on the mla by the tube lens (object side of tl)
    dio = obj2tube - dot; %% image side of the objective
    dof = Camera.fobj * dio/(dio - Camera.fobj); %% object side of the objective -> dof is focused on the MLA
    if isnan(dof) 
        dof = Camera.fobj;
    end
%     M_mla = tube2mla/dof;
    M_mla = Camera.M; %% in 4f systems magnification does not change with depth ?
    tubeRad = (dio-obj2tube) * objRad/dio;
    
    % if tube2mla is known and mla2sensor has to be retrived s.t. F number matching condition is satisfied
    if Camera.mla2sensor == 0
        mla2sensor = tube2mla*Camera.lensPitch/2/tubeRad;
        Camera.mla2sensor = mla2sensor;
    end
    
elseif (Camera.plenoptic == 1)
    
    if Camera.tube2mla == 0
        Camera.tube2mla = Camera.ftl;
    end
    
    if Camera.mla2sensor == 0
        Camera.mla2sensor = Camera.fm;
    end
    
    tubeRad = objRad;
    dof = Camera.fobj; %% object side of the objective -> dof is focused on the mla
    M_mla = Camera.M; % magnification up to the mla position
end

%% Write extra computed parameters to Camera struct
uRad =  tubeRad*Camera.mla2sensor/Camera.tube2mla;
offsetFobj = dof - Camera.fobj;
Camera.objRad = objRad;
Camera.uRad = uRad;
Camera.tubeRad = tubeRad;
Camera.dof = dof;
Camera.offsetFobj = offsetFobj;
Camera.M_mla = M_mla;