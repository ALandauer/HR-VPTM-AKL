

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing using Lagrangian tracking
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
if ~exist('results','dir') 
   mkdir('results')
end
v = VideoWriter('results/video_3D_inc_cum.avi'); v.FrameRate = 5; open(v); figure,
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



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Modify codes below to plot interpolated displacements and strains on a uniform grid mesh
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
Plotdisp_show3(uvw_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor');

%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show3(F_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor',1,tstep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% save out the results
%%

results_file_names_Lag = fullfile('results',['results_3D_LagTotal_',file_names{1}(1:end-4),'.mat']);
if ~exist('results','dir') 
   mkdir('results')
end

resultDispInc = resultDisp;
resultDefGradInc = resultDefGrad;
save(results_file_names_Lag,'resultDispInc','resultDefGradInc','parCoordTrajMat','disp_A2BCum','beadParam_all',...
    'MPTPara','F_Grid_refB_Vector','uvw_Grid_refB_Vector','coordinatesFEM_refB', 'elementsFEM_refB');

disp('%%%%% Cumulated reuslts saved %%%%%'); fprintf('\n');

%%
%plot mean strain comps
disp('%%%%% Compute tracked cumulative displacements %%%%%'); fprintf('\n');
parCoordTrajMat = cell2mat(parCoordTraj);
clear disp_A2BCum* parCoordA parCoordB parCoordACum parCoordBCum parPosHist

[row1,col1] = find(isnan(parCoordTrajMat(1:length(file_names):end,1))==0);
trackParCum_ind = row1;
trackParCum_track_ratio = [];

deform_func = -[0.025:0.025:0.25]; 
% expt_strain22 = deform_func-1;
% expt_strain11 = 1./sqrt(deform_func)-1;
% expt_strain33 = 1./sqrt(deform_func)-1;

for ImgSeqNum = 2:length(file_names)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(file_names):end,1))==0);
    trackParCum_ind = intersect(row2,trackParCum_ind);
    trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) / size(parCoord_prev{1},1);
    
    parCoordA = parCoordTrajMat(1:length(file_names):end,1:3);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(file_names):end,1:3);
    
    parCoordACum{ImgSeqNum-1} = parCoordA(trackParCum_ind,:);
    parCoordBCum{ImgSeqNum-1} = parCoordB(trackParCum_ind,:);
    disp_A2BCum_exp{ImgSeqNum-1} = parCoordBCum{ImgSeqNum-1} - parCoordACum{ImgSeqNum-1};
    
    %     parNum = find(trackParCum_ind == 10);
    %     parPosHist(ImgSeqNum-1) = parCoordBCum{ImgSeqNum-1}(parNum,3)';
    
end

clear mean_strain_* std_strain_* RMSD_strain_* u_total F_total

for ii = 1:length(disp_A2BCum_exp)
    
    disp_A2BCum_vec = disp_A2BCum_exp{ii};
    parCoordACum_vec = parCoordACum{ii};
    parCoordBCum_vec = parCoordBCum{ii};
%     parCoordBCum_vec = parCoordBCum_GT{ii};
    
    sxyz = min([round(0.5*MPTPara.f_o_s),20]);%.*MPTPara.axesScale; % Step size for griddata
    smoothness = 0.05; % Smoothness for regularization; "smoothness=0" means no regularization
    
    [x_Grid_refB,y_Grid_refB,z_Grid_refB,u_Grid_refB] = ...
        funScatter2Grid3D(parCoordACum_vec(:,1),parCoordACum_vec(:,2),parCoordACum_vec(:,3),disp_A2BCum_vec(:,1),sxyz,smoothness);
    [~,~,~,v_Grid_refB] = funScatter2Grid3D(parCoordACum_vec(:,1),parCoordACum_vec(:,2),parCoordACum_vec(:,3),disp_A2BCum_vec(:,2),sxyz,smoothness);
    [~,~,~,w_Grid_refB] = funScatter2Grid3D(parCoordACum_vec(:,1),parCoordACum_vec(:,2),parCoordACum_vec(:,3),disp_A2BCum_vec(:,3),sxyz,smoothness);
    
    % Build a displacement vector
    uvw_Grid_refB_Vector = [u_Grid_refB(:),v_Grid_refB(:),w_Grid_refB(:)]';
    uvw_Grid_refB_Vector = uvw_Grid_refB_Vector(:);
    % Calculate deformation gradient
    D_Grid = funDerivativeOp3(size(x_Grid_refB,1),size(x_Grid_refB,2),size(x_Grid_refB,3), sxyz ); % Central finite difference operator
    F_Grid_refB_Vector = D_Grid*uvw_Grid_refB_Vector; % {F} - I ={D}{U}
    
    %exp strain from def grad
    F_total{ii}{1,1} = reshape(F_Grid_refB_Vector(1:9:end),size(x_Grid_refB));
    F_total{ii}{2,2} = reshape(F_Grid_refB_Vector(5:9:end),size(x_Grid_refB));
    F_total{ii}{3,3} = reshape(F_Grid_refB_Vector(9:9:end),size(x_Grid_refB));
    F_total{ii}{1,2} = 2*reshape(0.5*(F_Grid_refB_Vector(2:9:end)+F_Grid_refB_Vector(4:9:end)),size(x_Grid_refB));
    F_total{ii}{1,3} = 2*reshape(0.5*(F_Grid_refB_Vector(3:9:end)+F_Grid_refB_Vector(7:9:end)),size(x_Grid_refB));
    F_total{ii}{2,3} = 2*reshape(0.5*(F_Grid_refB_Vector(6:9:end)+F_Grid_refB_Vector(8:9:end)),size(x_Grid_refB));
    F_total{ii}{2,1} = F_total{ii}{1,2}; F_total{ii}{3,1} = F_total{ii}{1,3}; F_total{ii}{3,2} = F_total{ii}{2,3};
    
end
%%
bd_wd = 5;
for ii = 1:length(F_total)
    mean_strain_11(ii,1) = mean(F_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_22(ii,1) = mean(F_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_33(ii,1) = mean(F_total{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_12(ii,1) = mean(F_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_13(ii,1) = mean(F_total{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_23(ii,1) = mean(F_total{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    
    std_strain_11(ii,1) = std(F_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
    std_strain_22(ii,1) = std(F_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
    std_strain_33(ii,1) = std(F_total{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
    std_strain_12(ii,1) = std(F_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
    std_strain_13(ii,1) = std(F_total{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
    std_strain_23(ii,1) = std(F_total{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
end

%u_total = inc2cum(u,dm,m,'linear');


%
step_num = 1:length(mean_strain_11);
figure
shadedErrorBar(step_num,mean_strain_11,std_strain_11,'g-x',1)
hold on
shadedErrorBar(step_num,mean_strain_22,std_strain_22,'b-*',1)
shadedErrorBar(step_num,mean_strain_33,std_strain_33,'r-+',1)
shadedErrorBar(step_num,mean_strain_12,std_strain_12,'y-.o',1)
shadedErrorBar(step_num,-mean_strain_13,std_strain_13,'m-.s',1)
shadedErrorBar(step_num,mean_strain_23,std_strain_23,'k-.^',1)
xlabel('step num')
ylabel('Strain')

title('Shear; Stdev shaded region')
