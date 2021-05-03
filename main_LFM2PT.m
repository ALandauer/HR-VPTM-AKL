%main run script to run LFM2PT: coolect data, configure and run reconstruction, run Trial-MPT, and save output.


% This project is based on two underlaying packages: oLaF for
% reconstruction and Trial-MPT for particle tracking.
%
% See:  A Stefanoiu, J Page, P Symvoulidis, GG Westmeyer, and T Lasser,
%       "Artifact-free deconvolution in light field microscopy,"
%       Opt. Express 27, 31644-31666 (2019)
%
%       and
%
%       [Jin's Trial-MPT paper]
%

clear, close all

%% Add dependecy dirs to path
import2ws();


%% Specify data to import
% data_folder = ['.',filesep,'data',filesep]; %main data directory
% data_subfolder = ['dz_002_20210415',filesep]; %subfolder for a specific experiment
% fileNamePrefix = 'dz_*';

data_folder = ['.',filesep,'data',filesep]; %main data directory
data_subfolder = ['exp_shear_2_45_to120',filesep]; %subfolder for a specific experiment
fileNamePrefix = '10sps_*';


%% ======================= SET-UP SECTION ============================

% ----------------------- oLaF parameters -----------------------

%reconstruction depth range (im um)
depthRange = [-650, 400];
% axial slice step (in um)
depthStep = 5;

% choose lenslet spacing (in  pixels) to downsample the number of pixels between mlens for speed up
newSpacingPx = 15; % 'default' means no up/down-sampling (newSpacingPx = lensPitch/pixelPitch)
% choose super-resolution factor as a multiple of lenslet resolution (= 1 voxel/lenslet)
superResFactor = 'default'; % default means sensor resolution

% ----------------------- Trial-MPT parameters -----------------------

%problem dimensions setup
MPTPara.DIM = 3;
MPTPara.axesScale = [1.64,1.64,depthStep]; % unit: um/px
MPTPara.depthRange = depthRange; % unit: um/px
MPTPara.tstep = 1; % unit: us

MPTPara.mode = 'inc'; % {'inc': incremental mode; 'cum': cumulative mode}

% Bead detection method
BeadPara.detectionMethod = 3; % {1-TPT code; 2-regionprops; 3-LFM}
% Image binary mask file
im_roi_mask_file_path = '';

%%  ======================= oLaF SECTION ============================

% get lenslet images for recon
[LensletImageSeq, imageNames, WhiteImage, configFile] = LFM_selectImages(data_folder,data_subfolder,fileNamePrefix);

figure; imagesc(LensletImageSeq{1}); colormap inferno; title ('LF image #1'); drawnow

% Specific LFM configuration and camera parameters (um units)
Camera = LFM_setCameraParams(configFile, newSpacingPx);

% Compute LFPSF Patterns and other prerequisites: lenslets centers, resolution related
[LensletCenters, Resolution, LensletGridModel, NewLensletGridModel] = ...
    LFM_computeGeometryParameters(Camera, WhiteImage, depthRange, depthStep, superResFactor, 1);

[H, Ht] = LFM_computeLFMatrixOperators(Camera, Resolution, LensletCenters);

% Correct the input image
% obtain the transformation between grid models
FixAll = LFM_retrieveTransformation(LensletGridModel, NewLensletGridModel);

% apply the transformation to the lenslet and white images
for ii = 1:length(LensletImageSeq)
    [correctedLensletImageSeq{ii}, correctedWhiteImage] = LFM_applyTransformation(LensletImageSeq{ii}, WhiteImage, FixAll, LensletCenters, 1);
    correctedLensletImageSeq{ii}(correctedLensletImageSeq{ii} < mean(correctedLensletImageSeq{ii}(:))) = mean(correctedLensletImageSeq{ii}(:));
    correctedLensletImageSeq{ii} = mat2gray(single(correctedLensletImageSeq{ii}));
end

% Reconstruct
for ii = 1:length(LensletImageSeq)
    % precompute image/volume sizes
    imgSize = size(correctedLensletImageSeq{ii});
    imgSize = imgSize + (1-mod(imgSize,2)); % ensure odd size
    
    texSize = ceil(imgSize.*Resolution.texScaleFactor);
    texSize = texSize + (1-mod(texSize,2)); % ensure odd size
    
    % Setup function pointers
    if (strcmp(Camera.focus, 'single'))
        backwardFUN = @(projection) LFM_backwardProject(Ht, projection, LensletCenters, Resolution, texSize, Camera.range);
        forwardFUN = @(object) LFM_forwardProject(H, object, LensletCenters, Resolution, imgSize, Camera.range);
    elseif (strcmp(Camera.focus, 'multi'))
        backwardFUN = @(projection) LFM_backwardProjectMultiFocus(Ht, projection, LensletCenters, Resolution, texSize, Camera.range);
        forwardFUN = @(object) LFM_forwardProjectMultiFocus(H, object, LensletCenters, Resolution, imgSize, Camera.range);
    else
        error('Invalid micro-lens type.')
    end
    
    % build anti-aliasing filter kernels
    lanczosWindowSize = 2;
    widths = LFM_computeDepthAdaptiveWidth(Camera, Resolution);
    lanczos2FFT = LFM_buildAntiAliasingFilter([texSize, length(Resolution.depths)], widths, lanczosWindowSize);
    
    % Apply EMS deconvolution
    it = 3;
    
    % initialization
    initVolume = ones([texSize, length(Resolution.depths)]);
    LFimage = correctedLensletImageSeq{ii};
    
    % background correction
    onesvol = ones(size(initVolume));
    onesForward = forwardFUN(onesvol);
    onesBack = backwardFUN(onesForward);
    
    reconVolume = deconvEMS(forwardFUN, backwardFUN, LFimage, it, initVolume, 1, lanczos2FFT, onesForward, onesBack);
    
    % Display the first reconstruction
% %     if ii == 1
% %         figure;
% %         if(size(reconVolume, 3) > 1)
% %             imshow3D(reconVolume, [], 'inferno');
% %         else
% %             imagesc3D(reconVolume); colormap inferno
% %         end
% %         
% %         
% %         %verify with user if the recon is okay
% %         prompt = '\nContinue deconvolving the remainder of the sequence? Y/N [Y]: ';
% %         yn = input(prompt,'s');
% %         if isempty(yn)
% %             yn = 'Y';
% %         end
% %         if strcmpi(yn,'y')
% %             disp('Onward!')
% %         else
% %             error('Deconv not accepted...')
% %         end
% %     end
    
    filename_cur = [imageNames{ii}(1:end-4),'.mat'];
    
    LensletImage = LensletImageSeq{ii};
    %save out the recon'd image
    save(filename_cur,'reconVolume','LensletImage','filename_cur')
    
end


%%  =================== TRIAL-MPT SECTION ============================

%summarize setup
disp('*************************************************************');
disp('Starting Trial-MPT...');
disp(['Tracking mode: ',MPTPara.mode]);
disp('*************************************************************'); fprintf('\n');


% get the data folder and name from the recon save
cur_dir = dir();
cur_dir = cur_dir(1).folder;
try
    [fileFolder,fileName,ext] = fileparts(filename_cur);
catch
    tempf = dir([data_folder,data_subfolder,'*.mat']);
    [fileFolder,~,~] = fileparts([tempf(1).folder,filesep,tempf(1).name]);
end

%%%%% Trial-MPT path %%%%%
fileTrialMPTPath = cur_dir;

%%%%% Particle detection parameters %%%%%
% Bead parameter setup
BeadPara.deconvIter = 5;           % Threshold for detecting particles
BeadPara.thres = 0.5;           % Threshold for detecting particles
BeadPara.beadSize = 30;          % Estimated radius of a single particle
BeadPara.minSize = 1000;           % Minimum num of voxels of a single particle
BeadPara.maxSize = 10000;         % Maximum num of voxels of a single particle
BeadPara.winSize = [50,50,50];     % By default
BeadPara.randNoise = 1e-7;      % By default
BeadPara.numBeadsPSF = 2;
BeadPara.fileFolder = fileFolder; %folder for raw images
BeadPara.saveIntermediates = 1;   %1 = save and reuse intermediate raw images
BeadPara.PSF = [];                % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
BeadPara.distMissing = 70;        % Distance threshold to check whether particle has a match or not 
BeadPara.color = 'white';         % By default

% Trial-MPT tracking

%%%%% Trial-MPT Parameter %%%%%
MPTPara.f_o_s = 500;              % Size of search field: max(|u|,|v|,|w|)
MPTPara.edge_width = 20;          % Width (in px) of image boarder within which to discard particles (due to poor tracking, and z-recon)
MPTPara.n_neighborsMax = 16;     % Max # of neighboring particles
MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
MPTPara.gbSolver = 2;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
MPTPara.smoothness = 0.05;       % Coefficient of regularization
MPTPara.outlrThres = 3;          % Threshold for removing outliers in MPT
MPTPara.maxIterNum = 8;         % Max ADMM iteration number
MPTPara.iterStopThres = 1e-3;    % ADMM iteration stopping threshold

MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;
MPTPara.post_proc_type = 'eulerian'; %Whether to use a Lagrangian trajectory routine or Eulerian gridded routine for post-proc

%%%% Postprocessing: merge trajectory segments %%%%%
%FOR LAGRANGIAN POST-PROCESSING 
MPTPara.strain_f_o_s = Inf;       % Size of virtual strain gauge
MPTPara.strain_n_neighbors = 16; % # of neighboring particles used in strain gauge
distThres = 55; % distance threshold to connect split trajectory segments
extrapMethod = 'pchip';  % extrapolation scheme to connect split trajectory segments
                         % suggestion: 'nearest' for Brownian motion
minTrajSegLength = 10;    % the minimum length of trajectory segment that will be extrapolated
maxGapTrajSeqLength = 1; % the max frame# gap between connected trajectory segments

%%%% Postprocessing: gridding %%%%%
%FOR EULERIAN POST-PROCESSING 
grid_spacing = [30,30,20]; % grid spacing parameters


%%%%% Run Trial-MPT tracking %%%%%
run_Trial_MPT_3D_hardpar_inc;

