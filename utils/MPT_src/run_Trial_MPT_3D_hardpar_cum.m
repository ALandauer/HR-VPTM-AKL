% %%%%%%%%%%%%%%%%%% Trial-MPT (3D cumulative mode) %%%%%%%%%%%%%%%%%
% Main file of code "Topology-based rotation-invariant augmented Lagrangian
% multiple particle tracking (Trial-MPT)"
% ***********************************************
% Dimension: 3D
% Tracking mode: cumulative deformation
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
% Author: Jin Yang
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Date: 2020.12.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%% Load 3D volumetric images %%%%%
try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end % Open image folder

ImgSeqNum=1; [file_names,Img] = funReadImage3([data_folder,data_subfolder,fileNamePrefix,'.mat'],ImgSeqNum); % Load image

try if isempty(fileFolder)~=1, cd(fileTrialMPTPath); end; catch; end % Come back to the main path

MPTPara.xRange = [0,size(Img{1},1)-1]*MPTPara.axesScale(1);
MPTPara.yRange = [0,size(Img{1},2)-1]*MPTPara.axesScale(1);

%%%%% Update MPTPara %%%%%
MPTPara.gridxyzROIRange.gridx = [1,size(Img{1},1)];
MPTPara.gridxyzROIRange.gridy = [1,size(Img{1},2)];
MPTPara.gridxyzROIRange.gridz = [1,size(Img{1},3)];

% figure, imagesc3(Img{1}) % To see volumetric image
disp('%%%%%% Load reference image: Done! %%%%%%'); fprintf('\n');

%%%%% Load image mask file %%%%%
try load(im_roi_mask_file_path); catch; end
try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if it exists
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
% BeadPara.maxSize = 1000;        % Maximum radius of a single particle
% BeadPara.winSize = [5, 5, 5];   % By default
% BeadPara.dccd = [1,1,1];        % By default
% BeadPara.abc = [1,1,1];         % By default
% BeadPara.forloop = 1;           % By default
% BeadPara.randNoise = 1e-7;      % By default
% BeadPara.numBeadsPSF = 1        % Number of bead to select and average for PSF
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
    BeadPara.detectionMethod = 3;
end
%%%%% Method 1: TPT code, better for high-density seeding %%%%%
if BeadPara.detectionMethod == 1
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    x_px{1}{ImgSeqNum} = locateParticles(double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:))),beadParam_all{ImgSeqNum}); % Detect particles
    x_sub{1}{ImgSeqNum} = radialcenter3dvec(double(Img{ImgSeqNum}),x_px{1}{ImgSeqNum},beadParam_all{ImgSeqNum}); % Localize particles
    x_sub{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units
% ----------------------------
    
%%%%% Method 2: Modified TracTrac code, better for lower density, medium size beads %%%%%
elseif BeadPara.detectionMethod == 2
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    x_sub{1}{ImgSeqNum} = f_detect_particles3(double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:))),beadParam_all{ImgSeqNum});
    x_sub{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units
    
%%%%% Method 3: Deconv + Active contour code, better for large beads that need bespoke deconv %%%%%
elseif BeadPara.detectionMethod == 3
    
    %method specific beadPara entries
    BeadPara.deconvThresh = 0.05;
    BeadPara.deconvPrefilter = true; %true/false gaussian prefilter option
    BeadPara.deconvIter = 5;
    BeadPara.psfSize = [25,25]; %x,y size of bead-based psf
    BeadPara.winSize = [7, 7, 7];
    BeadPara.ratThresh = 0.20;
    BeadPara.circThresh = 1.0;
    BeadPara.smoothFac = 0.15;
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    
    vol_in = double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:)));
    
    %run preprocessing to get PSF and deconvolve
    [vol_in,beadParam_all{ImgSeqNum}] = funPreprocLocalizeAC(vol_in,beadParam_all{ImgSeqNum},ImgSeqNum);
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
parCoord_prev = cell(length(file_names),1);     parCoord_prev{1} = parCoordA;
uvw_B2A_prev = cell(length(file_names)-1,1);    track_A2B_prev = cell(length(file_names)-1,1);
resultDisp = cell(length(file_names)-1,1);      resultDefGrad = cell(length(file_names)-1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2 : length(file_names)  % "ImgSeqNum" is the frame index
    
    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);
    
    %%%%% Load image volumetric data %%%%%
    try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end % Open image folder
    tempvol = load(file_names{ImgSeqNum}); fieldName = fieldnames(tempvol);
    Img{2} = getfield(tempvol,fieldName{2}); clear tempvol; %#ok<GFLD>
    if iscell(Img{2}), Img{2}=Img{2}{1}; end
    try if isempty(fileFolder)~=1, cd(fileTrialMPTPath); end; catch; end % Come back to the main path
    
    %%%%% Trial_MPT_tracking %%%%%
    [parCoordB_temp,uvw_B2A_temp,resultDisp_temp,resultDefGrad_temp,track_A2B_temp,~] = fun_TrialMPT_3D_HardPar( ...
        ImgSeqNum,Img{2},BeadPara,beadParam_all,MPTPara,parCoordA,parCoord_prev(2:end),uvw_B2A_prev);
    
    %%%%% Store results %%%%%
    parCoord_prev{ImgSeqNum} = parCoordB_temp;
    uvw_B2A_prev{ImgSeqNum-1} = uvw_B2A_temp;  % cumulative displacement
    resultDisp{ImgSeqNum-1} = resultDisp_temp;
    resultDefGrad{ImgSeqNum-1} = resultDefGrad_temp;
    track_A2B_prev{ImgSeqNum-1} = track_A2B_temp;
    
end


%%%%% Cumulative tracking ratio %%%%%
disp('%%%%% Calculate cumulative tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(file_names)-1,1);
DefType = 'exp'; defList = [2:1:length(file_names)]';

for ImgSeqNum = 2 : length(file_names)
    track_A2B = track_A2B_prev{ImgSeqNum-1};
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{ImgSeqNum},1);
end

fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Tracking ratio');
axis([2,length(file_names),0,1]);

%%%%% Save results %%%%%
disp('%%%%%% Trial-MPT hard particle tracking: Done! %%%%%%'); fprintf('\n');
results_file_name = 'results_3D_hardpar.mat';
save(results_file_name,'parCoord_prev','uvw_B2A_prev','resultDisp','resultDefGrad','track_A2B_prev');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Visualize tracked cumulative displacement of each frame %%%%%
disp('%%%%% Plot tracked cumulative deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
%already handled in localization
%try axesScale = MPTPara.xstep; catch, axes_scale = [1,1,1]; end % unit: um/px
axes_scale = [1,1,1];
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us
% ImgSeqNum  % Frame #

%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_3D_cum.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = 2:length(file_names)
    
    % Displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uvw_B2A_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};
    
    %%%%% Plot displacements %%%%%
    clf, plotCone3(parCoordB(:,1)*axes_scale(1),parCoordB(:,2)*axes_scale(2),parCoordB(:,3)*axes_scale, ...
        disp_A2B_parCoordB(:,1)*axes_scale(1),disp_A2B_parCoordB(:,2)*axes_scale(2),disp_A2B_parCoordB(:,3)*axes_scale(3));
    set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
    title(['Tracked cumulative displacement (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([MPTPara.xRange(1), MPTPara.xRange(2), ...
        MPTPara.yRange(1), MPTPara.yRange(2), ...
        MPTPara.depthRange(1), MPTPara.depthRange(2)]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Compute trajectory %%%%%

%%%%% Initialization %%%%%
resultDispCurr = resultDisp{1};
parCoordA = resultDispCurr.parCoordA;
parCoordATraj = cell(size(parCoordA,1),1);

%%%%% Compute and collect all trajectory segments %%%%%
for parInd = 1:size(parCoordA,1)
    
    for ImgSeqNum = 2:(size(resultDisp,1)+1)
        
        resultDispCurr = resultDisp{ImgSeqNum-1};
        parCoordB = resultDispCurr.parCoordB;
        track_A2B = resultDispCurr.track_A2B;
        
        if track_A2B(parInd) > 0
            parCoordATraj{parInd}(ImgSeqNum-1,1:3) = parCoordB(track_A2B(parInd),1:3);
        else
            parCoordATraj{parInd}(ImgSeqNum-1,1:3) = [nan,nan,nan];
        end
    end
    
end


%%%%% Plot tracked trajectories %%%%%
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');
figure,
for parInd = 1:size(parCoordA,1)
    try
        wayPoints = parCoordATraj{parInd};
        if (size(resultDisp,1)+1)<4
            hold on; line(wayPoints(isnan(wayPoints(:,1))<1,1),wayPoints(isnan(wayPoints(:,1))<1,2),wayPoints(isnan(wayPoints(:,1))<1,3)); view(3); % straight lines
        else
            hold on; fnplt(cscvn(wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1);
        end
        % if sum(1-isnan(wayPoints(:,1)))>1 % Don't show if there is only one point on the trajectory
        %    hold on; plot3(wayPoints(:,1),wayPoints(:,2),wayPoints(:,3),'.','markersize',8);
        % end
    catch
    end
end

set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
title('Tracked particle trajectory','fontweight','normal');
xlabel('x'); ylabel('y'); zlabel('z');
axis([MPTPara.xRange(1), MPTPara.xRange(2), ...
        MPTPara.yRange(1), MPTPara.yRange(2), ...
        MPTPara.depthRange(1), MPTPara.depthRange(2)]);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Modify codes below to plot interpolated displacements and strains on a uniform grid mesh');
pause;

ImgSeqNum = 2; % Frame #


%%%%% Previously tracked displacement field %%%%%
resultDispCurr = resultDisp{ImgSeqNum-1};
resultDefGradCurr = resultDefGrad{ImgSeqNum-1};
disp_A2B_parCoordB = resultDispCurr.disp_A2B_parCoordB;
parCoordB = resultDispCurr.parCoordB;

% % %%%%% remove rigid body translations %%%%%
% % disp_A2B_parCoordB(:,1) = disp_A2B_parCoordB(:,1) - median(disp_A2B_parCoordB(:,1));
% % disp_A2B_parCoordB(:,2) = disp_A2B_parCoordB(:,2) - median(disp_A2B_parCoordB(:,2));
% % disp_A2B_parCoordB(:,3) = disp_A2B_parCoordB(:,3) - median(disp_A2B_parCoordB(:,3));


%%%%% Interpolate scatterred data to gridded data %%%%%
sxyz = min([round(0.5*MPTPara.f_o_s),20]).*MPTpara.axesScale; % Step size for griddata
smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization

[x_Grid_refB,y_Grid_refB,z_Grid_refB,u_Grid_refB]=funScatter2Grid3D(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),disp_A2B_parCoordB(:,1),sxyz,smoothness);
[~,~,~,v_Grid_refB]=funScatter2Grid3D(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),disp_A2B_parCoordB(:,2),sxyz,smoothness);
[~,~,~,w_Grid_refB]=funScatter2Grid3D(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),disp_A2B_parCoordB(:,3),sxyz,smoothness);

% Build a displacement vector
uvw_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:),w_Grid_refB(:)]'; uvw_Grid_refB_Vector=uvw_Grid_refB_Vector(:);

% Calculate deformation gradient
D_Grid = funDerivativeOp3(size(x_Grid_refB,1),size(x_Grid_refB,2),size(x_Grid_refB,3),sxyz); % Central finite difference operator
F_Grid_refB_Vector=D_Grid*uvw_Grid_refB_Vector; % {F}={D}{U}


%%%%% Cone plot grid data: displecement %%%%%
figure, plotCone3(x_Grid_refB*axes_scale(1),y_Grid_refB*axes_scale(2),z_Grid_refB*axes_scale(3),u_Grid_refB*axes_scale(1),v_Grid_refB*axes_scale(2),w_Grid_refB*axes_scale(3));
set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
title('Tracked cumulative displacement','fontweight','normal');
axis([MPTPara.xRange(1), MPTPara.xRange(2), ...
        MPTPara.yRange(1), MPTPara.yRange(2), ...
        MPTPara.depthRange(1), MPTPara.depthRange(2)]);

%%%%% Generate an FE-mesh %%%%%
[coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp3(x_Grid_refB*axes_scale(1),y_Grid_refB*axes_scale(2),z_Grid_refB*axes_scale(3));

%%%%% Cone plot grid data: displacement %%%%%
Plotdisp_show3(uvw_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor');

%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show3(F_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor',1,tstep);













