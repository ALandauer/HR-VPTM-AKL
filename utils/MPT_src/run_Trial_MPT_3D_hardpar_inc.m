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

ImgSeqNum=1; [file_names,Img] = funReadImage3([data_folder,data_subfolder,fileNamePrefix,'.mat'],ImgSeqNum); % Load image

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
        ImgSeqNum,Img{2},BeadPara,beadParam_all,MPTPara,parCoord_prev{ImgSeqNum-1},parCoord_prev(2:end),uvw_B2A_prev);
    
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
results_file_names = fullfile('results',['results_3D_',file_names{1}(1:end-4),'.mat']);
save(results_file_names,'parCoord_prev','uvw_B2A_prev','track_A2B_prev','track_B2A_prev','resultDisp','resultDefGrad','beadParam_all');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Visualize tracked incremental displacement of each frame %%%%%
disp('%%%%% Plot tracked incremental deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
%already handled in localization
%try axes_scale = MPTPara.xstep; catch, axes_scale = [1,1,1]; end % unit: um/px
axes_scale = [1,1,1];
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us

%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_3D_inc.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = 2:length(file_names) % ImgSeqNum: Frame #
    
    % Displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uvw_B2A_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};
    
    %%%%% Plot displacements %%%%%
    clf, plotCone3(parCoordB(:,1)*axes_scale(1), parCoordB(:,2)*axes_scale(3), parCoordB(:,3)*axes_scale(3), ...
        disp_A2B_parCoordB(:,1)*axes_scale(1)/tstep, disp_A2B_parCoordB(:,2)*axes_scale(2)/tstep, disp_A2B_parCoordB(:,3)*axes_scale(3)/tstep );
    set(gca,'fontsize',18); box on; axis equal; view(3);
    title(['Tracked velocity (#',num2str(ImgSeqNum),')'],'fontweight','normal');
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
disp('%%%%% Merge tracked trajectory segments %%%%%'); fprintf('\n');

%%%%% Initialization %%%%%
%already handled in localization
%try axes_scale = MPTPara.xstep; catch, axes_scale = [1,1,1]; end % unit: um/px
axes_scale = [1,1,1];

try tstep = MPTPara.tstep; catch tstep = 1; end % time gap between consecutive frames
trajInd = 0; % index of trajectory segments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Compute and collect all trajectory segments %%%%%
for tempk = 1 : size(parCoord_prev,1)  % Find trajectories passing through particles in frame #tempk
    
    parCoordCurr = parCoord_prev{tempk}; % Find trajectories passing parCoordCurr
    clear parCoordTrajCurr; 
    parCoordTrajCurr = cell(size(parCoordCurr,1),1); % Initialize a cell structure to store #tempk trajectories
    
    % Add particles in frame #tempk to "parCoordTrajCurr"
    for tempj = 1:length(parCoordCurr)
        parCoordTrajCurr{tempj}(tempk,1:3) = parCoordCurr(tempj,1:3);
    end
    
    for ImgSeqNum = (1+tempk) : size(parCoord_prev,1) % For each later tracked incremental deformation
        
        parCoordB = parCoord_prev{ImgSeqNum}; % Current particle coordinates in frame #ImgSeqNum
        parCoordB_Ind = []; parCoordCurr_Ind = []; % Initialize particle index
        track_B2A_Curr = track_B2A_prev{ImgSeqNum-1}; % Current tracking in frame #ImgSeqNum
        for tempi = 1:length(track_B2A_Curr)
            try
                % Back propogation
                if ImgSeqNum > 2
                    for tempj = (ImgSeqNum-2) : -1 : tempk
                        track_B2A_Last = track_B2A_prev{tempj};
                        track_B2A_Curr(tempi) = track_B2A_Last(track_B2A_Curr(tempi));
                    end
                end
                % Check whether the trajectory has already been plotted previously
                if tempk>1
                    track_B2A_SecLast = track_B2A_prev{tempk-1};
                    track_B2A_SecLast_temp = track_B2A_SecLast(track_B2A_Curr(tempi)); % if this is not 0, means it's already been plotted
                else % Trajectories from first two frames (tempk=1) will be added
                    track_B2A_SecLast_temp = 0;
                end
                % Assign index values
                if (track_B2A_Curr(tempi)>0) && (track_B2A_SecLast_temp==0)
                    parCoordB_Ind = [parCoordB_Ind; tempi];
                    parCoordCurr_Ind = [parCoordCurr_Ind; track_B2A_Curr(tempi)];
                end
            catch
            end
        end
        
        for tempj = 1:length(parCoordCurr_Ind) % Add found trajectories to cell structure "parCoordTraj"
            parCoordTrajCurr{parCoordCurr_Ind(tempj)}(ImgSeqNum,1:3) = parCoordB(parCoordB_Ind(tempj),1:3);
        end
    end
    
    for parInd = 1:size(parCoordTrajCurr,1)
        wayPoints = parCoordTrajCurr{parInd};
        if ~isempty(wayPoints)
            wayPoints(wayPoints(:,1)==0,:) = wayPoints(wayPoints(:,1)==0,:)*nan;
            wayPoints = [wayPoints; nan(size(parCoord_prev,1)-size(wayPoints,1),3)];
            trajInd = trajInd + 1;
            parCoordTraj{trajInd} = wayPoints;
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Merge trajectory segments %%%%%
% Find the starting point and length of each trajectory segment
% parCoordTrajPara has size [xxx, 2], and each row is
% parCoordTrajPara(xxx, 1:2) = [first non-NAN position, trajectory segment length]
parCoordTrajPara = []; parCoordTraj = parCoordTraj(:);
for tempi = 1 : size(parCoordTraj,1)
    parCoordTrajCurr = parCoordTraj{tempi};
    parCoordTrajPara(tempi,1) = find(isnan(parCoordTrajCurr(:,1))==0, 1, 'first');
    parCoordTrajPara(tempi,2) = sum(1-isnan(parCoordTrajCurr(:,1)));
end

%%%%% Try to merge trajectory segments by extrapolating the particle position %%%%%
for tempMergeTime = 1:4 % Try to merge four times
    
    for tempm = 0:maxGapTrajSeqLength % tempm is the # of missing particles between trajectory segments
        
        for tempk = (size(parCoord_prev,1)-1) : -1 : minTrajSegLength % For trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]
            
            [row,col] = find( parCoordTrajPara(:,2)==tempk ); % Find trajectory segments whose "length==tempk"
            [row1,~] = find( parCoordTrajPara(:,2)<size(parCoord_prev,1)+1-tempk & parCoordTrajPara(:,2)>0 ); % Find trajectory setments with length requirement
            
            for tempi = 1:length(row) % For each trajectory segment whose "length==tempk"
                
                parCoordTrajMat = cell2mat( parCoordTraj );
                parCoordTrajCurr = parCoordTraj{row(tempi)}; % For each trajectory segment whose "length==tempk"
                parCoordTrajCurrPred_x = fillmissing(parCoordTrajCurr(:,1),extrapMethod); % fill missing using 'pchip' method
                parCoordTrajCurrPred_y = fillmissing(parCoordTrajCurr(:,2),extrapMethod); % fill missing using 'pchip' method
                parCoordTrajCurrPred_z = fillmissing(parCoordTrajCurr(:,3),extrapMethod); % fill missing using 'pchip' method
                % figure, plot3(parCoordTrajCurrPred_x,parCoordTrajCurrPred_y,parCoordTrajCurrPred_z,'-o');
                
                
                %%%%% Find all probable trajectory segments in the positive direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && (sum(parCoordTrajPara(row(tempi),1:2))+tempm>0)
                    
                    [row2,col2] = find( parCoordTrajPara(:,1) == sum(parCoordTrajPara(row(tempi),1:2))+tempm ); % starting point requirement
                    row3 = intersect(row1,row2);
                    temp1 = size(parCoord_prev,1)*(row3-1) +(sum(parCoordTrajPara(row(tempi),1:2)))+tempm; temp1=temp1'; temp1=temp1(:); % Index in parCoordTrajMat
                    temp2 = parCoordTrajMat(temp1, 1:3);
                    
                    % hold on; plot3(temp2(:,1),temp2(:,2),temp2(:,2),'.'); pause;
                    
                    temp3 = sqrt((temp2(:,1)-parCoordTrajCurrPred_x( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 + ...
                        (temp2(:,2)-parCoordTrajCurrPred_y( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 + ...
                        (temp2(:,3)-parCoordTrajCurrPred_z( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 );
                    
                    [temp3min,temp3minind] = min(temp3);
                    if temp3min < distThres % Find the continuing trajectory segment %JY!!!! threshold distance 3
                        
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 ) = ...
                            parCoordTraj{row3(temp3minind)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(temp3minind)} = parCoordTraj{row3(temp3minind)}*nan;
                        
                        % Update varaible "parCoordTrajPara" for trajectory segment length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(temp3minind),2) + tempm;
                        parCoordTrajPara(row3(temp3minind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        temp_z = fillmissing(temp(:,3),extrapMethod); % fill missing
                        
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3) = [temp_x, temp_y, temp_z];
                        
                    end
                    
                end
                
                
                %%%%% Find all probable trajectory segments in the negative direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && ((parCoordTrajPara(row(tempi),1)-1)-tempm>0)
                    
                    [row2,col2] = find( sum(parCoordTrajPara,2) == parCoordTrajPara(row(tempi),1)-tempm ); % ending point requirement
                    row3 = intersect(row1,row2);
                    temp1 = size(parCoord_prev,1)*(row3-1) + (parCoordTrajPara(row(tempi),1)-1)-tempm; temp1=temp1'; temp1=temp1(:); % Index in parCoordTrajMat
                    temp2 = parCoordTrajMat(temp1, 1:3);
                    
                    % hold on; plot(tempyy(:,1),tempyy(:,2),'.');
                    temp3 = sqrt( ( temp2(:,1)-parCoordTrajCurrPred_x( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 + ...
                        ( temp2(:,2)-parCoordTrajCurrPred_y( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 + ...
                        ( temp2(:,3)-parCoordTrajCurrPred_z( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2  );
                    
                    [temp3min,temp3minind] = min(temp3);
                    if temp3min < distThres % Find the continuing trajectory segment
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 ) = ...
                            parCoordTraj{row3(temp3minind)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(temp3minind)} = parCoordTraj{row3(temp3minind)}*nan;
                        % Update varaible "parCoordTrajPara" for both trajectory segment starting point and its length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(temp3minind),2) + tempm;
                        parCoordTrajPara(row(tempi),1) = parCoordTrajPara(row3(temp3minind),1);
                        parCoordTrajPara(row3(temp3minind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        temp_z = fillmissing(temp(:,3),extrapMethod); % fill missing
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3) = [temp_x, temp_y, temp_z];
                        
                    end
                    
                end
                
            end % End of "for each trajectory segment whose "length==tempk" "
            
        end % End of "for trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]"
        
    end % End of tempm
    
end % End of "for tempMergeTime = 1:5 % Try to merge four times"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot tracked trajectories %%%%%
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');

figure,
% CN = [36 63 146; 20 76 156; 28 100 175; 9 115 186; 23 128 196; 33 148 210;
%     25 158 218; 19 172 226; 25 186 229; 69 196 221; 118 205 214; 157 216 215;
%     169 220 217; 193 229 224; 216 237 233; 234 246 245]/255;
% CN = CN(end-size(parCoord_prev,1):end, 1:3); % Colormap I used in the manuscript


for tempi = 1:size(parCoordTrajPara,1)
    
    wayPoints = parCoordTraj{tempi};
    
    if (length(wayPoints(isnan(wayPoints(:,1))<1,1))+1)<4
        hold on; line(axes_scale(1)*wayPoints(isnan(wayPoints(:,1))<1,1), ...
            axes_scale(2)*wayPoints(isnan(wayPoints(:,1))<1,2), ...
            axes_scale(3)*wayPoints(isnan(wayPoints(:,1))<1,3), 'linewidth', 1); view(2); % straight lines
    else
        hold on; fnplt(cscvn([axes_scale.*wayPoints(isnan(wayPoints(:,1))<1,:)]'),'',1);
    end
    
    %%%%% Codes to plot trajectories with frame-dependent color %%%%%
    % if sum(1-isnan(wayPoints(:,1)))>1  % Don't show if there is only one point on the trajectory
    %     hold on; plot3(xstep*wayPoints(:,1),xstep*wayPoints(:,2),xstep*wayPoints(:,3),'r.','markersize',5);
    % end
    %
    % for tempj = 1:size(parCoord_prev,1)-1
    %     hold on; line(xstep*[wayPoints(tempj,1),wayPoints(tempj+1,1)], ...
    %         xstep*[wayPoints(tempj,2),wayPoints(tempj+1,2)], ...
    %         xstep*[wayPoints(tempj,3),wayPoints(tempj+1,3)], 'linewidth',1.2, 'color', CN(tempj,:) );
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
title('Tracked particle trajectory','fontweight','normal');
xlabel(''); ylabel(''); zlabel('z');
axis([MPTPara.xRange(1), MPTPara.xRange(2), ...
    MPTPara.yRange(1), MPTPara.yRange(2), ...
    MPTPara.depthRange(1), MPTPara.depthRange(2)]);



%%
%%%%% Compute cumulative tracking ratio from emerged trajectories %%%%%
disp('%%%%% Plot tracked cumulative displacements %%%%%'); fprintf('\n');
parCoordTrajMat = cell2mat( parCoordTraj );

[row1,col1] = find(isnan(parCoordTrajMat(1:length(file_names):end,1))==0);
trackParCum_ind = row1;
trackParCum_track_ratio = [];

for ImgSeqNum = 2:length(file_names)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(file_names):end,1))==0);
    trackParCum_ind = intersect(row2,trackParCum_ind);
    trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) / size(parCoord_prev{1},1);
end

defList = [2:1:length(file_names)];
fig=figure; ax=axes; hold on; plot(defList,trackParCum_track_ratio,'bs--','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Image #'); ylabel('Tracking ratio');
axis([2,length(file_names),0,1]);

%%%%% Plot tracked cumulative displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_3D_inc_cum.avi'); v.FrameRate = 5; open(v); figure,
for ImgSeqNum = 2:length(file_names)
    
    parCoordA = parCoordTrajMat(1:length(file_names):end,1:3);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(file_names):end,1:3);
    parCoordACum = parCoordA(trackParCum_ind,:);
    parCoordBCum = parCoordB(trackParCum_ind,:);
    disp_A2BCum = parCoordBCum - parCoordACum;
    
    % ----- Cone plot grid data: displecement -----
    clf; plotCone3(axes_scale(1)*parCoordBCum(:,1),axes_scale(2)*parCoordBCum(:,2),axes_scale(3)*parCoordBCum(:,3), ...
        axes_scale(1)*disp_A2BCum(:,1),axes_scale(2)*disp_A2BCum(:,2),axes_scale(3)*disp_A2BCum(:,3));
    set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
    title(['Tracked cumulative disp (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([MPTPara.xRange(1), MPTPara.xRange(2), ...
        MPTPara.yRange(1), MPTPara.yRange(2), ...
        MPTPara.depthRange(1), MPTPara.depthRange(2)]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);

%% %%% Compute cumulative displacement at each step %%%%%
disp('%%%%% Compute tracked cumulative displacements %%%%%'); fprintf('\n');
parCoordTrajMat = cell2mat(parCoordTraj);
clear disp_A2BCum RMSD_* parCoordA parCoordB parCoordACum parCoordBCum parPosHist

[row1,col1] = find(isnan(parCoordTrajMat(1:length(file_names):end,1))==0);
trackParCum_ind = row1;
trackParCum_track_ratio = [];

for ImgSeqNum = 2:length(file_names)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(file_names):end,1))==0);
    trackParCum_ind = intersect(row2,trackParCum_ind);
    trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) / size(parCoord_prev{1},1);
    
    parCoordA = parCoordTrajMat(1:length(file_names):end,1:3);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(file_names):end,1:3);
    parCoordACum{ImgSeqNum-1} = parCoordA(trackParCum_ind,:);
    parCoordBCum{ImgSeqNum-1} = parCoordB(trackParCum_ind,:);
    disp_A2BCum{ImgSeqNum-1} = parCoordBCum{ImgSeqNum-1} - parCoordACum{ImgSeqNum-1};
    
%     parNum = find(trackParCum_ind == 10);
%     parPosHist(ImgSeqNum-1) = parCoordBCum{ImgSeqNum-1}(parNum,3)'; 
    
end

% % % figure
% % % subplot(1,2,1)
% % % plot([-196:10:200],-(parPosHist))
% % % axis image
% % % subplot(1,2,2)
% % % plot([-196:10:200],-(parPosHist+[-196:10:200]))
% % % axis image

mean_cum_disp = cellfun(@(x) mean(x,1),disp_A2BCum,'UniformOutput',false);
mean_cum_disp = reshape(cell2mat(mean_cum_disp),3,[])';

std_cum_disp = cellfun(@(x) std(x,[],1),disp_A2BCum,'UniformOutput',false);
std_cum_disp = reshape(cell2mat(std_cum_disp),3,[])';


for ii = 1:length(disp_A2BCum)
    disp_meas_y_ = disp_A2BCum{ii}(:,1);
    disp_meas_y = disp_A2BCum{ii}(:,1);
    disp_meas_x = disp_A2BCum{ii}(:,2);
    disp_meas_z = disp_A2BCum{ii}(:,3);
    
%     disp_meas_y(abs(disp_meas_y_) > abs(mean(disp_meas_y_)+3*std(disp_meas_y_))) = [];
%     disp_meas_x(abs(disp_meas_y_) > abs(mean(disp_meas_y_)+3*std(disp_meas_y_))) = [];
%     disp_meas_z(abs(disp_meas_y_) > abs(mean(disp_meas_y_)+3*std(disp_meas_y_))) = [];
    
    N = length(disp_meas_z);
    
    disp_imps_x = zeros(N,1);
    disp_imps_y = zeros(N,1);
    
    disp_imps_z = -ii*20*ones(N,1);
    
    RMSD_y(ii,1) = sqrt(sum((disp_meas_y - disp_imps_y).^2)/N);
    RMSD_x(ii,1) = sqrt(sum((disp_meas_x - disp_imps_x).^2)/N);
    RMSD_z(ii,1) = sqrt(sum((disp_meas_z - disp_imps_z).^2)/N);
end

x_lbl = 'Experimental imposed disp z';
% imps_disp = 20*[1:length(mean_cum_disp)]';
imps_disp = -20*[1:length(mean_cum_disp)]';
% imps_disp = [0.022,0.025,0.028,0.033,0.040,0.050,0.066,0.100,0.200];
figure
subplot(1,3,1)
shadedErrorBar(imps_disp,mean_cum_disp(:,2),RMSD_x)
% shadedErrorBar(imps_disp,mean_cum_disp(:,2),std_cum_disp(:,2))
% xlabel('Noise level')
xlabel(x_lbl)
ylabel('Measured displacement in x, um')
% axis image

subplot(1,3,2)
shadedErrorBar(imps_disp,mean_cum_disp(:,1),RMSD_y)
% shadedErrorBar(imps_disp,mean_cum_disp(:,1),std_cum_disp(:,1))
xlabel(x_lbl)
% xlabel('Noise level')
ylabel('Measured displacement in y, um')
% axis image

subplot(1,3,3)
shadedErrorBar(imps_disp,mean_cum_disp(:,3),RMSD_z)
% shadedErrorBar(imps_disp,mean_cum_disp(:,3),std_cum_disp(:,3))
xlabel(x_lbl)
ylabel('Measured displacement in z, um')
title('RMSD shaded region')
% axis image

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Modify codes below to plot interpolated displacements and strains on a uniform grid mesh');
% pause;
clear disp_A2BCum parCoordA parCoordB parCoordACum parCoordBCum trackParCum_ind

ImgSeqNum = 2; % Frame #

parCoordTrajMat = cell2mat( parCoordTraj );

[row1,col1] = find(isnan(parCoordTrajMat(1:length(file_names):end,1)) == 0);
trackParCum_ind{1} = row1;

for ImgSeqNum_ = 2:length(file_names)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum_:length(file_names):end,1)) == 0);
    trackParCum_ind{ImgSeqNum_} = intersect(row2,trackParCum_ind{ImgSeqNum_-1});
end

%%%%% Previously tracked displacement field %%%%%
parCoordA = parCoordTrajMat(1:length(file_names):end,1:3);
parCoordB = parCoordTrajMat(ImgSeqNum:length(file_names):end,1:3);
parCoordACum = parCoordA(trackParCum_ind{ImgSeqNum},1:3);
parCoordBCum = parCoordB(trackParCum_ind{ImgSeqNum},1:3);
disp_A2BCum = parCoordBCum - parCoordACum;

%%%%% Interpolate scatterred data to gridded data %%%%%
sxyz = min([round(0.5*MPTPara.f_o_s),10]).*MPTPara.axesScale; % Step size for griddata
smoothness = 0.1;%1e-3; % Smoothness for regularization; "smoothness=0" means no regularization

[x_Grid_refB,y_Grid_refB,z_Grid_refB,u_Grid_refB]=funScatter2Grid3D(parCoordACum(:,1),parCoordACum(:,2),parCoordACum(:,3),disp_A2BCum(:,1),sxyz,smoothness);
[~,~,~,v_Grid_refB]=funScatter2Grid3D(parCoordACum(:,1),parCoordACum(:,2),parCoordACum(:,3),disp_A2BCum(:,2),sxyz,smoothness);
[~,~,~,w_Grid_refB]=funScatter2Grid3D(parCoordACum(:,1),parCoordACum(:,2),parCoordACum(:,3),disp_A2BCum(:,3),sxyz,smoothness);

% Apply ROI image mask
% [u_Grid_refB, v_Grid_refB] = funRmROIOutside(x_Grid_refB,y_Grid_refB,MPTPara.ImgRefMask,u_Grid_refB,v_Grid_refB);

% Build a displacement vector
uvw_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:),w_Grid_refB(:)]'; uvw_Grid_refB_Vector=uvw_Grid_refB_Vector(:);
% Calculate deformation gradient
D_Grid = funDerivativeOp3(size(x_Grid_refB,1),size(x_Grid_refB,2),size(x_Grid_refB,3), sxyz ); % Central finite difference operator
F_Grid_refB_Vector = D_Grid*uvw_Grid_refB_Vector; % {F}={D}{U}


%%%%% Cone plot grid data: displecement %%%%%
figure, plotCone3(axes_scale(1)*x_Grid_refB,axes_scale(2)*y_Grid_refB,axes_scale(3)*z_Grid_refB,u_Grid_refB*axes_scale(1),v_Grid_refB*axes_scale(2),w_Grid_refB*axes_scale(3));
set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
title('Tracked cumulative displacement','fontweight','normal');
axis([MPTPara.xRange(1), MPTPara.xRange(2), ...
    MPTPara.yRange(1), MPTPara.yRange(2), ...
    MPTPara.depthRange(1), MPTPara.depthRange(2)]);


%%%%% Generate an FE-mesh %%%%%
[coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp3(x_Grid_refB*axes_scale(1),y_Grid_refB*axes_scale(2),z_Grid_refB*axes_scale(3));

%%%%% Cone plot grid data: displacement %%%%%
%Plotdisp_show3(uvw_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor');

%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show3(F_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor',1,tstep);













