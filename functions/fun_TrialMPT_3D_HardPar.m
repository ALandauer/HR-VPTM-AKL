function [parCoordB,uvw_B2A_refB,resultDisp,resultDefGrad,track_A2B,track_B2A] = fun_TrialMPT_3D_HardPar(...
    ImgSeqNum,ImgDef,BeadPara,beadParam_all,MPTPara,parCoordA,parCoordB_prev,uvw_B2A_refB_prev)
 

warning('off');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: Particle detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Load current deformed frame %%%%%
currImg = ImgDef(MPTPara.gridxyzROIRange.gridx(1):MPTPara.gridxyzROIRange.gridx(2), ...
                 MPTPara.gridxyzROIRange.gridy(1):MPTPara.gridxyzROIRange.gridy(2), ...
                 MPTPara.gridxyzROIRange.gridz(1):MPTPara.gridxyzROIRange.gridz(2));

%%%%% If PSF is non-empty, perform deconvolution %%%%%
if ~isempty(BeadPara.PSF) && BeadPara.detectionMethod ~= 3
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp('----- Deconvolution frame #',num2str(ImgSeqNum),' ------');
end

%%%%% Pre-process bead image if bead color is "black" %%%%%
if strcmp(BeadPara.color,'black')
    disp('Not implement yet.')
else
    currImg2 = currImg;
end

%%%%% Several methods to detect particles %%%%%
try 
    BeadPara.detectionMethod = BeadPara.detectionMethod;
catch
    BeadPara.detectionMethod = 3;
end
%%%%% Method 1: TPT code %%%%%
if BeadPara.detectionMethod == 1
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    x{1}{ImgSeqNum} = locateParticles(double(currImg2)/max(double(currImg2(:))),beadParam_all{ImgSeqNum}); % Detect particles
    x{1}{ImgSeqNum} = radialcenter3dvec(double(currImg2),x{1}{ImgSeqNum},beadParam_all{ImgSeqNum}); % Localize particles
    x{1}{ImgSeqNum} = x{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units
% ----------------------------
%%%%% Method 2: Modified TracTrac code %%%%%
elseif BeadPara.detectionMethod == 2
    beadParam_all{ImgSeqNum} = funSetUpBeadParams(BeadPara);
    x{1}{ImgSeqNum} = f_detect_particles3(double(currImg2)/max(double(currImg2(:))),beadParam_all{ImgSeqNum});
    x{1}{ImgSeqNum} = x{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units

%%%%% Method 3: Deconv + Active contour code, better for large beads that need bespoke deconv %%%%%
elseif BeadPara.detectionMethod == 3
    
    %method specific beadPara entries
    BeadPara.deconvThresh = 0.05;
    BeadPara.deconvPrefilter = true; %true/false gaussian prefilter option
    BeadPara.deconvIter = 5;
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
    
    vol_in = double(currImg2)/max(double(currImg2(:)));
    
    %run preprocessing to get PSF and deconvolve
    [vol_in,beadParam_all{ImgSeqNum}] = funPreprocLocalizeAC(vol_in,beadParam_all{ImgSeqNum},ImgSeqNum);
    %find interger centriods
    [x_px{1}{ImgSeqNum},beadParam_all{ImgSeqNum}] = funLocateParticlesAC(vol_in,beadParam_all{ImgSeqNum},ImgSeqNum);
    %Use radial center-finding from TPT to get subpixel estimates based on the integer centroid locations
    x_sub{1}{ImgSeqNum} = radialcenter3dvec(double(currImg2),x_px{1}{ImgSeqNum},beadParam_all{ImgSeqNum});
    x_sub{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum}.*MPTPara.axesScale; %convert to um units
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Add MPTPara.gridxyzROIRange left-bottom corner point coordinates
x{1}{ImgSeqNum} = x_sub{1}{ImgSeqNum} + ...
    [MPTPara.gridxyzROIRange.gridx(1)*MPTPara.axesScale(1)+MPTPara.xRange(1)-1*MPTPara.axesScale(1), ...
    MPTPara.gridxyzROIRange.gridy(1)*MPTPara.axesScale(2)+MPTPara.yRange(1)-1*MPTPara.axesScale(2), ...
    MPTPara.gridxyzROIRange.gridz(1)*MPTPara.axesScale(3)+MPTPara.depthRange(1)-1*MPTPara.axesScale(3)];
parCoordB = x{1}{ImgSeqNum};

% Remove bad coordinates that are out of image ROI
%%%%% Remove parCoord outside the image area %%%%%
parCoordB( parCoordB(:,1) > MPTPara.xRange(2),:) = [];
parCoordB( parCoordB(:,2) > MPTPara.yRange(2),:) = [];
parCoordB( parCoordB(:,3) > MPTPara.depthRange(2),:) = [];
parCoordB( parCoordB(:,1) < MPTPara.xRange(1),:) = [];
parCoordB( parCoordB(:,2) < MPTPara.yRange(1),:) = [];
parCoordB( parCoordB(:,3) < MPTPara.depthRange(1),:) = [];

%%%%% Plot detected particles %%%%%
% close all;
% figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'ro'); 
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); 
% title('Detected particles in ref image','fontweight','normal');
% figure, plot3(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),'ro'); 
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); 
% title('Detected particles in def image','fontweight','normal');


%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp(['Detected particle # in def image: ',num2str(size(parCoordB,1))]);
% disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2: Particle Linking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Track A & B neighbor-topology match ======
%
%  Directly run this section if coordinates of detected particles are known:
%  Coordinates of particles in reference image:   parCoordA
%  Coordinates of particles in deformed image:    parCoordB
%
%    \  |  /                  \  |  /
%     \ | /                    \ | /
%   --- A ---       v.s.     --- B ---
%      / \                      / \
%     /   \                    /   \
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
[matches_A2B,u_B2A_curr_refB,track_A2B] = f_track_trial_match3D( parCoordA, parCoordB, ...
   'f_o_s',MPTPara.f_o_s, 'n_neighborsMax', MPTPara.n_neighborsMax, 'n_neighborsMin', MPTPara.n_neighborsMin, ...
   'gbSolver', MPTPara.gbSolver, 'smoothness', MPTPara.smoothness, ...
   'outlrThres', MPTPara.outlrThres, 'maxIterNum', MPTPara.maxIterNum, ...
   'iterStopThres', MPTPara.iterStopThres, 'usePrevResults', MPTPara.usePrevResults, ...
   'strain_n_neighbors',MPTPara.strain_n_neighbors, 'strain_f_o_s',MPTPara.strain_f_o_s, ...
   'gridxyzROIRange',MPTPara.gridxyzROIRange, 'parCoordB_prev',parCoordB_prev, ...
   'uvw_B2A_prev',uvw_B2A_refB_prev, 'ImgSeqNum',ImgSeqNum,'ImgSize',size(currImg), ...
   'BeadParaDistMissing',BeadPara.distMissing);


%%%%% Compute track_B2A %%%%%
matches_A2B = matches_A2B(matches_A2B(:,2)>0,:); % Build matches_A2B_temp
track_B2A = zeros(size(parCoordB,1), 1);
for tempi = 1:size(matches_A2B)
    track_B2A(matches_A2B(tempi,2)) = matches_A2B(tempi,1);
end  

 
%% %%%%% Plotting %%%%%%
% Compute displacement from tracked particles on deformed frame
disp_A2B_parCoordB = -u_B2A_curr_refB;
figure, plotCone3(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3), ...
 disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2),disp_A2B_parCoordB(:,3));
set(gca,'fontsize',18); box on; axis equal; axis tight; view(3);  
title('Tracked displacements','fontweight','normal');
xlabel(''); ylabel(''); cb = colorbar; set(cb,'fontsize',18);



%% %%%%% Compute F = def grad = grad_u = (grad_x X - I) %%%%%

% strain_f_o_s: size of virtual strain gauge
% strain_n_neighbors: # of neighboring particles used in strain gauge

%%%%% Strain based on Moving Least Square Fitting in deformed configuration %%%%%
[XYZ_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad3(disp_A2B_parCoordB, parCoordB, MPTPara.strain_f_o_s, MPTPara.strain_n_neighbors);

%%%%% Strain based on Moving Least Square Fitting in reference configuration %%%%%
[XYZ_refA,U_A2B_refA,F_A2B_refA] = funCompDefGrad3(disp_A2B_parCoordB, parCoordB-disp_A2B_parCoordB, MPTPara.f_o_s, MPTPara.strain_n_neighbors);

 

%% %%%%% Store results %%%%%
uvw_B2A_refB = u_B2A_curr_refB;

resultDisp = struct('parCoordA',parCoordA,'parCoordB',parCoordB,'track_A2B', ...
                         track_A2B,'disp_A2B_parCoordB',disp_A2B_parCoordB);
                     
resultDefGrad = struct('XY_refA',XYZ_refA,'U_A2B_refA',U_A2B_refA,'F_A2B_refA',F_A2B_refA, ...
                        'XY_refB',XYZ_refB,'U_B2A_refB',U_B2A_refB,'F_B2A_refB',F_B2A_refB);






