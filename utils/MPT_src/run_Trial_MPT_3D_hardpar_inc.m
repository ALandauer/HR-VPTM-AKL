% %%%%%%%%%%%%%%%%%% Trial-MPT (3D incremental mode) %%%%%%%%%%%%%%%%%
% Main file of code "Topology-based rotation-invariant augmented Lagrangian
% multiple particle tracking (Trial-MPT)"
% ***********************************************
% Dimension: 3D
% Tracking mode: incremental deformation
% Particle type: hard (particle shape is rigid)
% -----------------------------------------------
%
% -----------------------------------------------
% References
% [1] M Patel, SE Leggett, AK Landauer, IY Wong, C Franck. Rapid,
%     topology-based particle tracking for high-resolution measurements of
%     large complex 3D motion fields. Scientific Reports. 8:5581 (2018).
% [2] J Yang, L Hazlett, AK Landauer, C Franck. Augmented Lagrangian
%     Digital Volume Correlation (ALDVC). Experimental Mechanics (2020).
% [3] T Janke, R Schwarze, K Bauer. Part2Track: A MATLAB package for double
%     frame and time resolved Particle Tracking Velocimetry. 11, 100413, SoftwareX (2020).
% [4] J Heyman. TracTrac: a fast multi-object tracking algorithm for motion
%     estimation. Computers & Geosciences, vol 128, 11-18 (2019).
% [5] https://www.mathworks.com/matlabcentral/fileexchange/77347-gridded-interpolation-and-gradients-of-3d-scattered-data
% [6] https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% -----------------------------------------------
% Author: Jin Yang (jyang526@wisc.edu) and Alex Landauer (landauer@brown.edu)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%% Load 3D volumetric images %%%%%
try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end % Open image folder

ImgSeqNum=1; 
cur_image = [data_folder,data_subfolder,fileNamePrefix,'.mat'];
[file_names,Img] = funReadImage3(cur_image,ImgSeqNum); % Load image

try if isempty(fileFolder)~=1, cd(fileTrialMPTPath); end; catch; end % Come back to the main path

MPTPara.xRange = [MPTPara.edge_width,size(Img{1},1)-1-MPTPara.edge_width]*MPTPara.axesScale(1);
MPTPara.yRange = [MPTPara.edge_width,size(Img{1},2)-1-MPTPara.edge_width]*MPTPara.axesScale(2);

%%%%% Update MPTPara %%%%%
MPTPara.gridxyzROIRange.gridx = [1,size(Img{1},1)];
MPTPara.gridxyzROIRange.gridy = [1,size(Img{1},2)];
MPTPara.gridxyzROIRange.gridz = [1,size(Img{1},3)];

% figure, imagesc3(Img{1}) % To see volumetric image
disp('%%%%%% Load reference image: Done! %%%%%%'); fprintf('\n');

%%%%% Load image mask file %%%%%
try load(im_roi_mask_file_path); catch; end
try
    MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
catch
    disp('No mask, using whole image...')
    MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
end
disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');


%% ====== Detect particles ======
%%%%% Particle detection parameters %%%%%
%%%%% Bead Parameter %%%%%
% BeadPara.thres = 0.4;           % Threshold for detecting particles
% BeadPara.beadSize = 0;          % Estimated radius of a single particle
% BeadPara.minSize = 2;           % Minimum radius of a single particle
% BeadPara.maxSize = 1000;        % Maximum radius of a single particle
% BeadPara.winSize = [5, 5, 5];   % By default
% BeadPara.dccd = [1,1,1];        % By default
% BeadPara.abc = [1,1,1];         % By default
% BeadPara.forloop = 1;           % By default
% BeadPara.randNoise = 1e-7;      % By default
% BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
% BeadPara.distMissing = 5;       % Distance threshold to check whether particle has a match or not
% BeadPara.color = 'white';       % By default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgSeqNum = 1; % First reference image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Several methods to detect particles %%%%%
try
    BeadPara.detectionMethod = BeadPara.detectionMethod;
catch
    BeadPara.detectionMethod = 2;
end
%%%%% Method 1: TPT code %%%%%
if BeadPara.detectionMethod == 1
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    x_px{1}{ImgSeqNum} = locateParticles(double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:))),beadParam_all{ImgSeqNum}); % Detect particles
    x_sub{1}{ImgSeqNum} = radialcenter3dvec(double(Img{ImgSeqNum}),x_px{1}{ImgSeqNum},beadParam_all{ImgSeqNum}); % Localize particles
    x_sub{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units
    % ----------------------------
    %%%%% Method 2: Modified TracTrac code %%%%%
elseif BeadPara.detectionMethod == 2
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    x_sub{1}{ImgSeqNum} = f_detect_particles3(double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:))),beadParam_all{ImgSeqNum});
    x_sub{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units
    
    %%%%% Method 3: Deconv + Active contour code, better for large beads that need bespoke deconv %%%%%
elseif BeadPara.detectionMethod == 3
    
    %method specific beadPara entries
    BeadPara.deconvThresh = 0.05;
    BeadPara.deconvPrefilter = true; %true/false gaussian prefilter option
    BeadPara.psfSize = [50,50]; %x,y size of bead-based psf
    BeadPara.winSize = [51, 51, 51];
    BeadPara.ratThresh = 0.20;
    BeadPara.circThresh = 1.2;
    BeadPara.smoothFac = 0.15;
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    if ImgSeqNum > 1
        beadParam_all{ImgSeqNum}.minSize = beadParam_all{1}.minSize;
        beadParam_all{ImgSeqNum}.maxSize = beadParam_all{1}.maxSize;
        beadParam_all{ImgSeqNum}.thres = beadParam_all{1}.thres;
    end
    
    vol_in = double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:)));
    
    %run preprocessing to get PSF and deconvolve
    [vol_in,beadParam_all{ImgSeqNum}] = funPreprocLocalizeAC(vol_in,beadParam_all{ImgSeqNum},file_names,ImgSeqNum);
    %find interger centriods
    [x_px{1}{ImgSeqNum},beadParam_all{ImgSeqNum}] = funLocateParticlesAC(vol_in,beadParam_all{ImgSeqNum},ImgSeqNum);
    %Use radial center-finding from TPT to get subpixel estimates based on the integer centroid locations
    x_sub{1}{ImgSeqNum} = radialcenter3dvec(double(Img{ImgSeqNum}),x_px{1}{ImgSeqNum},beadParam_all{ImgSeqNum});
    x_sub{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Store particle positions as "parCoordA" %%%%%
x{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum} + ...
    [MPTPara.gridxyzROIRange.gridx(1)*MPTPara.axesScale(1)+MPTPara.xRange(1)-1*MPTPara.axesScale(1), ...
    MPTPara.gridxyzROIRange.gridy(1)*MPTPara.axesScale(2)+MPTPara.yRange(1)-1*MPTPara.axesScale(2), ...
    MPTPara.gridxyzROIRange.gridz(1)*MPTPara.axesScale(3)+MPTPara.depthRange(1)-1*MPTPara.axesScale(3)];
parCoordA = x{1}{ImgSeqNum};

%%%%% Remove parCoord outside the image area %%%%%
parCoordA( parCoordA(:,1) > MPTPara.xRange(2),:) = [];
parCoordA( parCoordA(:,2) > MPTPara.yRange(2),:) = [];
parCoordA( parCoordA(:,3) > MPTPara.depthRange(2),:) = [];
parCoordA( parCoordA(:,1) < MPTPara.xRange(1),:) = [];
parCoordA( parCoordA(:,2) < MPTPara.yRange(1),:) = [];
parCoordA( parCoordA(:,3) < MPTPara.depthRange(1),:) = [];

%%%%% Plot %%%%%
figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'bo');
view(3); box on; axis equal; axis tight; set(gca,'fontsize',18);
title('Detected particles in ref image','fontweight','normal');

%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');

%% %%%%% Initialization %%%%%
%%%%% MPT Parameter %%%%%
% MPTPara.f_o_s = 60;              % Size of search field: max(|u|,|v|,|w|)
% MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
% MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
% MPTPara.gbSolver = 2;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
% MPTPara.smoothness = 1e-1;       % Coefficient of regularization
% MPTPara.outlrThres = 5;          % Threshold for removing outliers in TPT
% MPTPara.maxIterNum = 20;         % Max ADMM iteration number
% MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
% MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
% MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge
% MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;


%%%%%% To store results %%%%%
parCoord_prev = cell(length(file_names)-1,1);  parCoord_prev{1} = parCoordA;
track_A2B_prev = cell(length(file_names)-1,1); track_B2A_prev = cell(length(file_names)-1,1);
uvw_B2A_prev = cell(length(file_names)-1,1);
resultDisp = cell(length(file_names)-1,1);
resultDefGrad = cell(length(file_names)-1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2 : length(file_names)  % "ImgSeqNum" is the frame index
    
    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);
    
    %%%%% Load image volumetric data %%%%%
    try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end % Open image folder
    tempvol = load(file_names{ImgSeqNum}); fieldName = 'reconVolume';
    Img{2} = getfield(tempvol,fieldName); clear tempvol; %#ok<GFLD>
    if iscell(Img{2}), Img{2}=Img{2}{1}; end
    try if isempty(fileFolder)~=1, cd(fileTrialMPTPath); end; catch; end % Come back to the main path
    
    %%%%% Trial_MPT_tracking %%%%%
    [parCoordB_temp,uvw_B2A_temp,resultDisp{ImgSeqNum-1},resultDefGrad{ImgSeqNum-1},track_A2B_temp,track_B2A_temp,beadParam_all] = fun_TrialMPT_3D_HardPar( ...
        ImgSeqNum,Img{2},file_names,BeadPara,beadParam_all,MPTPara,parCoord_prev{ImgSeqNum-1},parCoord_prev(2:end),uvw_B2A_prev);
    
    %%%%% Store results %%%%%
    parCoord_prev{ImgSeqNum} = parCoordB_temp;
    uvw_B2A_prev{ImgSeqNum-1} = uvw_B2A_temp; % incremental displacement
    track_A2B_prev{ImgSeqNum-1} = track_A2B_temp;
    track_B2A_prev{ImgSeqNum-1} = track_B2A_temp;
    
end


%%%%% Incremental tracking ratio %%%%%
disp('%%%%% Calculate incremental tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(file_names)-1,1);
defList = [2:1:length(file_names)]';

for ImgSeqNum = 2 : length(file_names)
    track_A2B = track_A2B_prev{ImgSeqNum-1};
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{ImgSeqNum},1);
end

fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Tracking ratio');
axis([2,length(file_names),0,1]);

%%%%% Save results %%%%%
disp('%%%%%% Trial-MPT 3D hard particle tracking: Done! %%%%%%'); fprintf('\n');
results_file_names_inc = fullfile('results',['results_3D_inc_',file_names{1}(1:end-4),'.mat']);
if ~exist('results','dir') 
   mkdir('results')
end

save(results_file_names_inc,'resultDisp','resultDefGrad','beadParam_all','MPTPara');

if strcomi(MPTPara.post_proc_type,'lagrangian') 
    run_post_process_lagrangian
elseif strcmpi(MPTPara.post_proc_type,'eulerian')
    run_post_process_eulerian
else
    disp('No post processing used')
end






