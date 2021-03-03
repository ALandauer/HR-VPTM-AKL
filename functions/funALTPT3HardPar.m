function [parCoordB_prevCurr,uvw_B2A_prevCurr,resultDispCurr,resultDefGradCurr,track_A2B] = funALTPT3HardPar(...
    ImgSeqNum,ImgDef,BeadPara,MPTPara,parCoordA,parCoordB_prev,uvw_B2A_prev)



warning('off');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: Particle detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Load current image frame %%%%%
currImg = ImgDef(MPTPara.gridxyzROIRange.gridx(1):MPTPara.gridxyzROIRange.gridx(2), ...
                 MPTPara.gridxyzROIRange.gridy(1):MPTPara.gridxyzROIRange.gridy(2), ...
                 MPTPara.gridxyzROIRange.gridz(1):MPTPara.gridxyzROIRange.gridz(2));

%%%%% If PSF is non-empty, perform deconvolution %%%%%
if ~isempty(BeadPara.PSF)
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp('----- Deconvolution frame #',num2str(ImgSeqNum),' ------');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pre-process bead image if bead color is black %%%%%
if strcmp(BeadPara.color,'black')
    disp('Not implement yet.')
    % ImgGauss = imgaussfilt(imgaussfilt(currImg,1),1);
    % ImgGauss(ImgGauss > BeadPara.thres*max(double(currImg(:)))) = 0;
    % bw = imbinarize(uint8(ImgGauss),'adaptive','ForegroundPolarity','dark','Sensitivity',0.8); % figure, imshow(bws2);
    % bws2 = bwareaopen(bw,round(pi*BeadPara.minSize^2)); % remove all object containing fewer than BeadPara.minSize
    % removeobjradius = BeadPara.minSize; % fill a gap in the pen's cap
    % se = strel('disk',removeobjradius);
    % bws2 = imclose(bws2,se);
    % currImg2 = double(bws2); % figure, imshow(uint8(ImgGauss));
else
    currImg2 = currImg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Several methods to detect particles %%%%%
%%%%% Method 1: TPT code %%%%%
x{1}{ImgSeqNum} = locateParticles(double(currImg2)/max(double(currImg2(:))),BeadPara); % Detect particles
x{1}{ImgSeqNum} = radialcenter3dvec(double(currImg2),x{1}{ImgSeqNum},BeadPara); % Localize particles
% ----------------------------
%%%%% Method 2: Modified TracTrac code %%%%%
% x{1}{ImgSeqNum} = f_detect_particles3(double(currImg2)/max(double(currImg2(:))),BeadPara);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% Add back MPTPara.gridxyzROIRange
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyzROIRange.gridx(1)-1, MPTPara.gridxyzROIRange.gridy(1)-1, MPTPara.gridxyzROIRange.gridz(1)-1];
parCoordB = x{1}{ImgSeqNum};

% Remove bad coordinates outside image ROI
for tempi=1:3, parCoordB( parCoordB(:,tempi)>size(ImgDef,tempi), : ) = []; end
for tempi=1:3, parCoordB( parCoordB(:,tempi)<1, : ) = []; end


% %%%%% Plot detected particles %%%%%
% close all;
% figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'ro'); 
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); 
% title('Detected particles in ref image','fontweight','normal');
% figure, plot3(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),'ro'); 
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); 
% title('Detected particles in def image','fontweight','normal');

%%%%% Report detected beads # %%%%%
% disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp(['Detected particle # in def image: ',num2str(size(parCoordB,1))]);
% disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2: Particle Matching
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

%%%%% MPTPara Initialization %%%%%
try f_o_s = MPTPara.f_o_s;                   catch f_o_s = 20;          end   % window size used in global step solver
try n_neighborsMax = MPTPara.n_neighborsMax; catch n_neighborsMax = 25; end   % Max # of neighboring particles
try n_neighborsMin = MPTPara.n_neighborsMin; catch n_neighborsMin = 1;  end   % Min # of neighboring particles
try Subpb2Solver = MPTPara.Subpb2Solver;     catch Subpb2Solver = 1;    end   % ALTPT global step solver: 0-moving least square fitting; 1-global regularization
try smoothness = MPTPara.smoothness;         catch smoothness = 1e-2;   end
try outlrThres = MPTPara.outlrThres;         catch outlrThres = 5;      end   % Threshold for removing outliers in TPT
try maxIterNum = MPTPara.maxIterNum;         catch maxIterNum = 20;     end   % Max iteration number
try iterThres = MPTPara.iterThres;           catch iterThres = 1e-2;    end
try gauss_interp = MPTPara.gauss_interp;     catch gauss_interp = 1;    end
try usePrevResults = MPTPara.usePrevResults; catch usePrevResults = 0;  end
try strain_n_neighbors = MPTPara.strain_n_neighbors; catch strain_n_neighbors = 20; end
try strain_f_o_s = MPTPara.strain_f_o_s;     catch strain_f_o_s = 60;   end
try iterStopThres = MPTPara.iterStopThres;   catch iterStopThres = 1e-2; end

sxyz = min([round(0.5*f_o_s),20])*[1,1,1];
 

iterNum = 0; tic; % Initialize iteration number and start timing
parCoordBCurr = parCoordB; u_B2A_curr_refB = 0*parCoordBCurr; % BCurr is the warped ImgB
AParNotMissing = [1:size(parCoordA,1)]'; BCurrParNotMissing = [1:size(parCoordBCurr,1)]'; % Remove missing particles in both ref and def images
MatchRatioEqualsOneTime = 0;


%%%%% Use last step's result for a later frame %%%%%
if ImgSeqNum>2 && usePrevResults
    disp('Use previous results as an initial estimate');
    [tempu,tempv,tempw] = funInitGuess3(parCoordB_prev,uvw_B2A_prev,parCoordBCurr,ImgSeqNum);
    u_B2A_curr_refB = u_B2A_curr_refB + [tempu,tempv,tempw];
    parCoordBCurr = parCoordBCurr + [tempu,tempv,tempw];
end

%%%%% Plotting %%%%%
% ===== Original image pair =====
% figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'bo'); view(2);
% hold on, plot3(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),'k+'); view(2);
% box on; axis equal; axis tight; set(gca,'fontsize',18); box on;
% title('Original particle position','fontweight','normal');
% lgd = legend('Particle A','Particle B','location','southoutside','fontsize',12);
% 
% % ===== Warped image pair =====
% figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'bo'); view(2);
% hold on, plot3(parCoordBCurr(:,1),parCoordBCurr(:,2),parCoordBCurr(:,3),'k+'); view(2);
% box on; axis equal; axis tight; set(gca,'fontsize',18); box on;
% title('Warped particle position','fontweight','normal');
% lgd = legend('Particle A','Particle BCurr','location','southoutside','fontsize',12);
% 
% pause;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while iterNum < maxIterNum
    
    close all;
    n_neighbors = round(n_neighborsMin + exp(-0.5*iterNum)*(n_neighborsMax-n_neighborsMin)); % # of neighboring particles: n_neighbors0 (topology) --> 1 (nearest neighboring search)
    iterNum = iterNum+1; disp(['------ Iter #',num2str(iterNum),' ------']); % current iteration number
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Local step: Neighbor topology match %%%%%
    matches_A2B = [];  
    while isempty(matches_A2B) && n_neighborsMax < length(AParNotMissing)
        
        if n_neighbors > 3
            matches_A2B = f_track_neightopo_match3( parCoordA(AParNotMissing,:), parCoordBCurr(BCurrParNotMissing,:), f_o_s, n_neighbors );
            % matches_A2B = f_track_hist_match3( parCoordA(AParNotMissing,:), parCoordBCurr(BCurrParNotMissing,:), f_o_s, n_neighbors, gauss_interp );
        else
            matches_A2B = f_track_nearest_neighbour3( parCoordA(AParNotMissing,:), parCoordBCurr(BCurrParNotMissing,:), f_o_s );
        end
        if isempty(matches_A2B) == 1
            n_neighborsMax = round(n_neighborsMax + 5);
        else
            matches_A2B = [AParNotMissing(matches_A2B(:,1)), BCurrParNotMissing(matches_A2B(:,2))];
            [track_A2B, u_A2B] = funCompDisp3(parCoordA, parCoordBCurr, matches_A2B, outlrThres);
            
            MatchRatio = size(matches_A2B,1) / length(AParNotMissing);
            disp( ['Tracking ratio: ', num2str( size(matches_A2B,1) ),'/', num2str(length(AParNotMissing)), ' = ', num2str(MatchRatio) ]);
            
        end
    end
    if isempty(matches_A2B) == 1, disp('Wrong!'); break; end
    
    
    % if iterNum>1
    %     plotCone2(parCoordB(:,1),parCoordB(:,2),u_B2A_curr(:,1),u_B2A_curr(:,2));
    %     set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  view(2);
    %     title('Iterative tracked disp field','fontweight','normal');
    % end
    % 
    % plotCone2(parCoordA(track_A2B>0,1),parCoordA(track_A2B>0,2),u_A2B(:,1),u_A2B(:,2)); view(2);
    % set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  view(2);
    % title('Iterative tracked disp field','fontweight','normal');
    % pause;
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Global step: Compute kinematically compatible displacement %%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Method I: local moving least squares %%%%%
    if Subpb2Solver==1
        
        [XYZ_B2A_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad3(-u_A2B, parCoordBCurr(track_A2B(track_A2B>0),:), f_o_s, n_neighbors);
        
        [row,~] = find(isnan(U_B2A_refB(1:3:end)) == 1);        % Find nans
        XYZ_B2A_refB(row,:) = [];                               % remove nans
        U_B2A_refB([3*row-2; 3*row-1; 3*row]) = [];             % remove nans
        F_B2A_refB([9*row-8; 9*row-7; 9*row-6; 9*row-5;  ...
            9*row-4; 9*row-3; 9*row-2; 9*row-1; 9*row]) = [];   % remove nans
        
        Fx = scatteredInterpolant(XYZ_B2A_refB,U_B2A_refB(1:3:end),'linear','linear');
        tempu = Fx(parCoordBCurr);
        Fy = scatteredInterpolant(XYZ_B2A_refB,U_B2A_refB(2:3:end),'linear','linear');
        tempv = Fy(parCoordBCurr);
        Fz = scatteredInterpolant(XYZ_B2A_refB,U_B2A_refB(3:3:end),'linear','linear');
        tempw = Fz(parCoordBCurr);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Method II: global regularization %%%%%
    elseif Subpb2Solver==2
        
        tempUVW_B2A = -u_A2B;
        tempXYZ_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        
        try smoothness = smoothness; % Smoothness for regularization; "smoothness=0" means no regularization
        catch smoothness = 1e-3; % By default
        end
        try
            try
                [xGrid_refB,yGrid_refB,zGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,smoothness);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,smoothness);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,smoothness);
            catch
                [xGrid_refB,yGrid_refB,zGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,0);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,0);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,0);
            end
        catch
            [xGrid_refB,yGrid_refB,zGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),[1,1,1],0);
            [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),[1,1,1],0);
            [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),[1,1,1],0);
        end
        Fx = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],uGrid_B2A_refB_iter(:),'linear','linear');
        tempu = Fx(parCoordBCurr);
        Fy = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],vGrid_B2A_refB_iter(:),'linear','linear');
        tempv = Fy(parCoordBCurr);
        Fz = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],wGrid_B2A_refB_iter(:),'linear','linear');
        tempw = Fz(parCoordBCurr);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Method III: augmented Lagrangian regularization %%%%%
    elseif Subpb2Solver==3
        
        tempUVW_B2A = -u_A2B;
        tempXYZ_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        
        try smoothness = smoothness; % Smoothness for regularization; "smoothness=0" means no regularization
        catch smoothness = 1e-3; % By default
        end
        if iterNum == 1
            [xGrid_refB,yGrid_refB,zGrid_refB] = ndgrid(MPTPara.gridxyzROIRange.gridx(1) : sxyz(1) : MPTPara.gridxyzROIRange.gridx(2), ...
                                                        MPTPara.gridxyzROIRange.gridy(1) : sxyz(2) : MPTPara.gridxyzROIRange.gridy(2), ...
                                                        MPTPara.gridxyzROIRange.gridz(1) : sxyz(3) : MPTPara.gridxyzROIRange.gridz(2) );
        end
        try 
            try
                [~,~,~,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,smoothness,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,smoothness,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,smoothness,xGrid_refB,yGrid_refB,zGrid_refB);
            catch
                [~,~,~,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,0,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,0,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,0,xGrid_refB,yGrid_refB,zGrid_refB);
            end
        catch
            [~,~,~,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),[1,1,1],0,xGrid_refB,yGrid_refB,zGrid_refB);
            [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),[1,1,1],0,xGrid_refB,yGrid_refB,zGrid_refB);
            [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),[1,1,1],0,xGrid_refB,yGrid_refB,zGrid_refB);
        end
        [M,N,L] = size(uGrid_B2A_refB_iter);
        DMat = funDerivativeOp3(M,N,L,mean(sxyz)*[1,1,1]);
        uVec = [uGrid_B2A_refB_iter(:),vGrid_B2A_refB_iter(:),wGrid_B2A_refB_iter(:)]'; uVec=uVec(:);
        
        if iterNum==1
            vdualVec=0*uVec; mu=1; alphaList=[1e-2,1e-1,1e0,1e1,1e2,1e3];
            
            % close all; [coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp(xGrid,yGrid);
            for tempi = 1:length(alphaList)
                alpha = alphaList(tempi);
                uhatVec = (alpha*(DMat')*DMat + speye(3*M*N*L))\((uVec-vdualVec));
                % Plotdisp_show(uhatVec,coordinatesFEM_refB,elementsFEM_refB);
                Err1(tempi) = sqrt((uhatVec-uVec+vdualVec)'*(uhatVec-uVec+vdualVec));
                Err2(tempi) = sqrt((DMat*uhatVec)'*(DMat*uhatVec));
            end
            
            ErrSum = Err1/max(Err1) + Err2/max(Err2);
            [~,indexOfalpha] = min(ErrSum);
            try % Tune the best beta by a quadratic polynomial fitting
                [fitobj] = fit(log10(alphaList(indexOfalpha-1:1:indexOfalpha+1))',ErrSum(indexOfalpha-1:1:indexOfalpha+1),'poly2');
                p = coeffvalues(fitobj); alpha_best = 10^(-p(2)/2/p(1));
            catch, alpha_best = alphaList(indexOfalpha);
            end
        end
        
        %%%%% Resolve global step with tuned best alpha %%%%%
        uhatVec = (alpha_best*(DMat')*DMat + speye(3*M*N*L))\((uVec-vdualVec));
        
        %%%%% Update dual variable
        vdualVec = vdualVec + uhatVec - uVec;
        
        %%%%% Interpolate to scatterred data points %%%%% 
        uGrid_B2A_refB_iter = reshape( uhatVec(1:3:end), size(uGrid_B2A_refB_iter) );
        vGrid_B2A_refB_iter = reshape( uhatVec(2:3:end), size(uGrid_B2A_refB_iter) );
        wGrid_B2A_refB_iter = reshape( uhatVec(3:3:end), size(uGrid_B2A_refB_iter) );
        
        Fx = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],uGrid_B2A_refB_iter(:),'linear','linear');
        tempu = Fx(parCoordBCurr);
        Fy = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],vGrid_B2A_refB_iter(:),'linear','linear');
        tempv = Fy(parCoordBCurr);
        Fz = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],wGrid_B2A_refB_iter(:),'linear','linear');
        tempw = Fz(parCoordBCurr);
        
        
    end % END of if Subpb2Solver==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Check convergence %%%%%
    u_B2ACurr_updateNorm = sqrt( ( norm(tempu)^2+norm(tempv)^2+norm(tempw)^2 )/length(tempu) );
    disp(['Disp update norm: ',num2str(u_B2ACurr_updateNorm)]);
    
    
    %%%%% If MatchRatio==1, only do additional five iterations %%%%%
    if MatchRatio == 1
        MatchRatioEqualsOneTime = MatchRatioEqualsOneTime+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if u_B2ACurr_updateNorm < sqrt(3)*iterStopThres || MatchRatioEqualsOneTime>5
        disp(['----- Converged! ------']); 
        break
        
    else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%% Warp back parCoordBCurr %%%%%
        u_B2A_curr_refB = u_B2A_curr_refB + [tempu,tempv,tempw];
        parCoordBCurr = parCoordBCurr + [tempu,tempv,tempw];
        
        % %%%%% Plot %%%%%
        % % ===== Original image pair =====
        % figure, plot(parCoordA(:,1),parCoordA(:,2),'bo'); view(2);
        % hold on; plot(parCoordB(:,1),parCoordB(:,2),'r+'); view(2);
        % box on; axis equal; axis tight; set(gca,'fontsize',18); box on;
        % title('Original particle position','fontweight','normal');
        % lgd = legend('Particle A','Particle B','location','southoutside','fontsize',12);
        % % ===== Matched image pair =====
        % figure, plot(parCoordA(matches_A2B(:,1),1),parCoordA(matches_A2B(:,1),2),'bo'); view(2);
        % hold on; plot(parCoordBCurr(matches_A2B(:,2),1),parCoordBCurr(matches_A2B(:,2),2),'r+'); view(2);
        % box on; axis equal; axis tight; set(gca,'fontsize',18); box on;
        % title('Current particle position','fontweight','normal');
        % lgd = legend('Matched Particle A','Matched Particle BCurr','location','southoutside','fontsize',12);
        % % ===== Warped image pair =====
        % figure, plot(parCoordA(:,1),parCoordA(:,2),'bo'); view(2);
        % hold on, plot(parCoordBCurr(:,1),parCoordBCurr(:,2),'k+'); view(2);
        % box on; axis equal; axis tight; set(gca,'fontsize',18); box on;
        % title('Warped particle position','fontweight','normal');
        % lgd = legend('Particle A','Particle BCurr','location','southoutside','fontsize',12);
        
        if n_neighbors < n_neighborsMax/2 
            %%%%% To detect missing particles in both reference and deformed images %%%%%
            neighborInd_BCurrA = knnsearch(parCoordBCurr(:,1:3),parCoordA(:,1:3),'K',1); % Find pts in BCurr near pts in ACurr
            dist_BCurrA = sqrt( sum((parCoordBCurr(neighborInd_BCurrA,1:3) -  parCoordA).^2,2) ); % figure, h = histogram(dist_BCurrA);
            % [AParNotMissing,~] = find(dist_BCurrA < max([2, max(abs(u_B2A_curr(:)))])); % Find particles not missing in Particle A
            [AParNotMissing,~] = find(dist_BCurrA < max([2, 0.5*BeadPara.distMissing])); % Find particles not missing in Particle A
            
            neighborInd_ABCurr = knnsearch(parCoordA(:,1:3),parCoordBCurr(:,1:3),'K',1); % Find pts in ACurr near pts in BCurr
            dist_ABCurr = sqrt( sum((parCoordA(neighborInd_ABCurr,1:3) -  parCoordBCurr).^2,2) );
            % [BCurrParNotMissing,~] = find(dist_ABCurr < max([2, max(abs(u_B2A_curr(:)))])); % Find particles not missing in Particle BCurr
            [BCurrParNotMissing,~] = find(dist_ABCurr < max([2, 0.5*BeadPara.distMissing]));  % Find particles not missing in Particle BCurr
        end
        
        %%%%% Improved warped image pair by removing missing particles %%%%%
%         figure; plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'kx'); view(2);
%         hold on, plot3(parCoordA(AParNotMissing,1),parCoordA(AParNotMissing,2),parCoordA(AParNotMissing,3),'ro'); view(2);
%         hold on, plot3(parCoordBCurr(:,1),parCoordBCurr(:,2),parCoordBCurr(:,3),'b+'); view(2);
%         hold on, plot3(parCoordBCurr(BCurrParNotMissing,1),parCoordBCurr(BCurrParNotMissing,2),parCoordBCurr(BCurrParNotMissing,3),'go'); view(2);
%         box on; axis equal; axis tight; set(gca,'fontsize',18); box on;
%         title('Warped particle position','fontweight','normal');
%         lgd = legend('Particle A','Partcile A not missing', ...
%             'Particle BCurr','Partcile BCurr not missing','location','southoutside','fontsize',12);
%         
%         pause;
        %%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%% Neighbor topology tracking: Done! %%%%%%'); fprintf('\n'); % timeCost = toc; toc


%% %%%%% Plotting %%%%%%
% Compute displacement from tracked particles on deformed frame
disp_A2B_parCoordB = -u_B2A_curr_refB;
figure,plotCone3(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2),disp_A2B_parCoordB(:,3));
set(gca,'fontsize',18); view(3); box on; axis equal; axis tight; view(2);
title('Tracked displacements','fontweight','normal');
xlabel(''); ylabel(''); title('');
cb = colorbar; set(cb,'fontsize',18);


% %%%%% Compute F = def grad = grad_u = (grad_x X - I) %%%%%
% strain_n_neighbors = 20; strain_f_o_s = 60;

% %%%%% Strain based on Moving Least Square Fitting in deformed configuration %%%%%
[XYZ_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad3(disp_A2B_parCoordB, parCoordB, strain_f_o_s, strain_n_neighbors);

% %%%%% Strain based on Moving Least Square Fitting in reference configuration %%%%%
[XYZ_refA,U_A2B_refA,F_A2B_refA] = funCompDefGrad3(disp_A2B_parCoordB, parCoordB-disp_A2B_parCoordB, f_o_s, strain_n_neighbors);




%% %%%%% Plot %%%%%
%     figure, imshow(imread(file_name{1,1}))
%     hold on; plot(parCoordA(:,1),parCoordA(:,2),'bo');
%     view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); title('Ref image','fontweight','normal');
%
%     figure, imshow(imread(file_name{1,1}))
%     hold on; plot( parCoordB(:,1)-disp_A2B(:,1), parCoordB(:,2)-disp_A2B(:,2), 'ro');
%     view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); title('Solved solution','fontweight','normal');

%%%%% Store results %%%%%
parCoordB_prevCurr = parCoordB;
uvw_B2A_prevCurr = u_B2A_curr_refB;

resultDispCurr = struct('parCoordA',parCoordA,'parCoordB',parCoordB,'track_A2B', ...
                         track_A2B,'disp_A2B_parCoordB',disp_A2B_parCoordB);
                     
resultDefGradCurr = struct('XYZ_refA',XYZ_refA,'U_A2B_refA',U_A2B_refA,'F_A2B_refA',F_A2B_refA, ...
                        'XYZ_refB',XYZ_refB,'U_B2A_refB',U_B2A_refB,'F_B2A_refB',F_B2A_refB);






