
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use a scattered2grid and incremental finite deformation multiplicative cumulation
%

% set up vars
smoothness = MPTPara.smoothness;
grid_spacing = [20,20,15]; 


%define data range, pad edges to accond for bead near edges moving outside
%the range of the particles in the first image
coords_cur = resultDisp{1}.parCoordA;
xList = MPTPara.xRange(1)-1.5*MPTPara.edge_width*MPTPara.axesScale(1):grid_spacing(1):MPTPara.xRange(2)+1.5*MPTPara.edge_width*MPTPara.axesScale(1);
yList = MPTPara.yRange(1)-1.5*MPTPara.edge_width*MPTPara.axesScale(2):grid_spacing(2):MPTPara.yRange(2)+1.5*MPTPara.edge_width*MPTPara.axesScale(2);
%z range is slightly different since it can be negative or both pos and neg
zRange_grid(1) = sign(min(coords_cur(:,3)))*(min(abs(coords_cur(:,3)))-2*grid_spacing(3));
zRange_grid(2) = sign(max(coords_cur(:,3)))*(max(abs(coords_cur(:,3)))+2*grid_spacing(3));
zList = min(zRange_grid):grid_spacing(3):max(zRange_grid);
[yGrid,xGrid,zGrid] = meshgrid(yList,xList,zList);

%Set z range without padding
MPTPara.zRange(1) = sign(min(coords_cur(:,3)))*(min(abs(coords_cur(:,3))));
MPTPara.zRange(2) = sign(max(coords_cur(:,3)))*(max(abs(coords_cur(:,3))));

disp('%%%%% Interpolating tracking results on ref grid %%%%%'); fprintf('\n');
%interpolate scattered data onto eluerian grid
u_inc = cell(1,length(resultDisp));
H_inc = cell(1,length(resultDefGrad));
for ii = 1:length(resultDisp)
    
    if ~mod(ii,10)
        disp(ii)
    end
    
    coords_cur_ = resultDisp{ii}.parCoordB;
    keep_coords = ~[coords_cur_(:,1) > MPTPara.xRange(2)|...
                   coords_cur_(:,2) > MPTPara.yRange(2)|...
                   coords_cur_(:,1) < MPTPara.xRange(1)|...
                   coords_cur_(:,2) < MPTPara.yRange(1)|...
                   abs(coords_cur_(:,3)) < min(abs(MPTPara.zRange))|...
                   abs(coords_cur_(:,3)) > max(abs(MPTPara.zRange))];
    coords_cur_disp = coords_cur_(keep_coords,:);
    
    disp_cur = resultDisp{ii}.disp_A2B_parCoordB(keep_coords,:);
    
    [x_grid_ref,y_grid_ref,z_grid_ref,u_inc{ii}{1}] = ...
                           funScatter2Grid3D(coords_cur_disp(:,1),coords_cur_disp(:,2),coords_cur_disp(:,3),disp_cur(:,1),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,u_inc{ii}{2}] = funScatter2Grid3D(coords_cur_disp(:,1),coords_cur_disp(:,2),coords_cur_disp(:,3),disp_cur(:,2),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,u_inc{ii}{3}] = funScatter2Grid3D(coords_cur_disp(:,1),coords_cur_disp(:,2),coords_cur_disp(:,3),disp_cur(:,3),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    
    coords_cur_ = resultDefGrad{ii}.XY_refA;
    keep_coords = ~[coords_cur_(:,1) > MPTPara.xRange(2)|...
                    coords_cur_(:,2) > MPTPara.yRange(2)|...
                    coords_cur_(:,1) < MPTPara.xRange(1)|...
                    coords_cur_(:,2) < MPTPara.yRange(1)|...
                    abs(coords_cur_(:,3)) < min(abs(MPTPara.zRange))|...
                    abs(coords_cur_(:,3)) > max(abs(MPTPara.zRange))];
    coords_cur_strain = coords_cur_(keep_coords,:);
                                         
    
    F_cur_vec = resultDefGrad{ii}.F_A2B_refA;
    F_cur_ = zeros(length(coords_cur_),6);
    F_cur_(:,1) = F_cur_vec(1:9:end);
    F_cur_(:,2) = F_cur_vec(5:9:end);
    F_cur_(:,3) = F_cur_vec(9:9:end);
    F_cur_(:,4) = F_cur_vec(4:9:end);
    F_cur_(:,5) = F_cur_vec(7:9:end);
    F_cur_(:,6) = F_cur_vec(8:9:end);
    F_cur = F_cur_(keep_coords,:)/2;
    
    [~,~,~,H_inc{ii}{1,1}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,1),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{2,2}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,2),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{3,3}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,3),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{1,2}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,4),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{1,3}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,5),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{2,3}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,6),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    H_inc{ii}{2,1} = H_inc{ii}{1,2}; H_inc{ii}{3,1} = H_inc{ii}{1,3}; H_inc{ii}{3,2} = H_inc{ii}{2,3};
end

%%
% compute the cumulative displacement field
disp('%%%%% Cumulating displacements %%%%%'); fprintf('\n');
msh{1} = xGrid;
msh{2} = yGrid;
msh{3} = zGrid;
u_total = inc2cum(u_inc,grid_spacing(1),msh,'linear');


%%
disp('%%%%% Cumulating deformation gradient %%%%%'); fprintf('\n');
% compute the cumulative deformation gradiant (F = e - I)
F_total = cell(length(H_inc)+1,1);
F_total{1} = cell(3);

for j = 1:3
    for k = 1:3
        if k == j
            F_total{1}{j,k} = ones(size(H_inc{1}{j,k}));
        else
            F_total{1}{j,k} = zeros(size(H_inc{1}{j,k}));
        end
    end
end

for ii = 1:length(H_inc)
    if ~mod(ii,10)
        disp(ii)
    end
    for loc = 1:numel(H_inc{ii}{1,1})
        H_inc_pt(loc,1) = H_inc{ii}{1,1}(loc);
        H_inc_pt(loc,2) = H_inc{ii}{2,2}(loc);
        H_inc_pt(loc,3) = H_inc{ii}{3,3}(loc);
        H_inc_pt(loc,4) = H_inc{ii}{1,2}(loc);
        H_inc_pt(loc,5) = H_inc{ii}{1,3}(loc);
        H_inc_pt(loc,6) = H_inc{ii}{2,3}(loc);
        
        %Convert to stretch
        Fhat_inc(1,1,loc) = H_inc_pt(loc,1)+1;
        Fhat_inc(2,2,loc) = H_inc_pt(loc,2)+1;
        Fhat_inc(3,3,loc) = H_inc_pt(loc,3)+1;
        
        Fhat_inc(1,2,loc) = H_inc_pt(loc,4);
        Fhat_inc(2,1,loc) = H_inc_pt(loc,4);
        Fhat_inc(1,3,loc) = H_inc_pt(loc,5);
        Fhat_inc(3,1,loc) = H_inc_pt(loc,5);
        Fhat_inc(2,3,loc) = H_inc_pt(loc,6);
        Fhat_inc(3,2,loc) = H_inc_pt(loc,6);
        
        [i,j,k] = ind2sub(size(xGrid),loc);
        if ii == 1
            %F_total = Finc for first step
            for m = 1:3
                for n = 1:3
                    F_total{ii+1}{m,n}(i,j,k) = Fhat_inc(m,n,loc);
                end
            end
        else
            %current F
            for m = 1:3
                for n = 1:3
                    Fcur(m,n,loc) = F_total{ii}{m,n}(loc);
                end
            end
            
            %F_total = F0*F1*F2...Fn
            F_tot_ = Fcur(:,:,loc)*Fhat_inc(:,:,loc);
            
            for m = 1:3
                for n = 1:3
                    F_total{ii+1}{m,n}(i,j,k) = F_tot_(m,n);
                end
            end
        end
    end
    
end

%% save the results
results_file_names_Eul = fullfile('results',['results_3D_EulTotal_',file_names{1}(1:end-4),'.mat']);
if ~exist('results','dir') 
   mkdir('results')
end

resultDispInc = resultDisp;
resultDefGradInc = resultDefGrad;
save(results_file_names_Eul,'resultDispInc','resultDefGradInc','u_total','F_total','beadParam_all','MPTPara');

disp('%%%%% Cumulated Eulerian reuslts saved %%%%%'); fprintf('\n');

%% plotting options
disp('%%%%% Optional plotting code - modify as needed %%%%%'); fprintf('\n');
% % % %%%%% Cone plot grid data: displecement %%%%%
% % % frame_num = 5;
% % % 
% % % figure, plotCone3(x_grid_ref,y_grid_ref,z_grid_ref,u_total{frame_num}{1},u_total{frame_num}{2},u_total{frame_num}{3});
% % % set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
% % % title('Tracked cumulative displacement','fontweight','normal');
% % % axis([MPTPara.xRange(1), MPTPara.xRange(2), ...
% % %     MPTPara.yRange(1), MPTPara.yRange(2), ...
% % %     MPTPara.depthRange(1), MPTPara.depthRange(2)]);
% % % 
% % % % uvw_Grid_ref_vector
% % % uvw_Grid_ref_vector=[u_total{frame_num}{1}(:),u_total{frame_num}{2}(:),u_total{frame_num}{3}(:)]'; 
% % % uvw_Grid_ref_vector=uvw_Grid_ref_vector(:);
% % % 
% % % % F_vector = [F11_pt1,F21_pt1,F31_pt1,F12_pt1,F22_pt1,F32_pt1,F13_pt1,F23_pt1,F33_pt1, ...
% % % %                   F11_pt2,F21_pt2,F31_pt2,F12_pt2,F22_pt2,F32_pt2,F13_pt2,F23_pt2,F33_pt2, ...
% % % %                   ... Other points ... ]'
% % % n_pts = numel(F_total{1}{1});
% % % F_grid_ref_vector = zeros(9*n_pts,1);
% % % F_grid_ref_vector(1:9:end) = F_total{frame_num}{1,1}(:);
% % % F_grid_ref_vector(2:9:end) = F_total{frame_num}{2,1}(:);
% % % F_grid_ref_vector(3:9:end) = F_total{frame_num}{3,1}(:);
% % % F_grid_ref_vector(4:9:end) = F_total{frame_num}{1,2}(:);
% % % F_grid_ref_vector(5:9:end) = F_total{frame_num}{2,2}(:);
% % % F_grid_ref_vector(6:9:end) = F_total{frame_num}{3,2}(:);
% % % F_grid_ref_vector(7:9:end) = F_total{frame_num}{1,3}(:);
% % % F_grid_ref_vector(8:9:end) = F_total{frame_num}{2,3}(:);
% % % F_grid_ref_vector(9:9:end) = F_total{frame_num}{3,3}(:);
% % % 
% % % %%%%% Generate an FE-mesh %%%%%
% % % [coordinatesFEM_ref,elementsFEM_ref] = funMeshSetUp3(x_grid_ref,y_grid_ref,z_grid_ref);
% % % 
% % % %%%%% Cone plot grid data: displacement %%%%%
% % % % Plotdisp_show3(uvw_Grid_ref_vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor');
% % % % 
% % % % %%%%% Cone plot grid data: infinitesimal strain %%%%%
% % % % Plotstrain_show3(F_grid_ref_vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor',1,tstep);



%%
%plot mean strain comps
bd_wd = 2;
clear mean_strain_* std_strain_*

[E_total,e_total] = calculateEij(F_total);

for ii = 1:length(E_total)-1

   
    
    mean_strain_11(ii,1) = mean(E_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_22(ii,1) = mean(E_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_33(ii,1) = mean(E_total{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_12(ii,1) = mean(E_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_13(ii,1) = mean(E_total{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    mean_strain_23(ii,1) = mean(E_total{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
    
    std_strain_11(ii,1) = std(E_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan')/sqrt(sum(track_A2B_prev{ii}>0));
    std_strain_22(ii,1) = std(E_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan')/sqrt(sum(track_A2B_prev{ii}>0));
    std_strain_33(ii,1) = std(E_total{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan')/sqrt(sum(track_A2B_prev{ii}>0));
    std_strain_12(ii,1) = std(E_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2))/2,[],'all','omitnan')/sqrt(sum(track_A2B_prev{ii}>0));
    std_strain_13(ii,1) = std(E_total{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2))/2,[],'all','omitnan')/sqrt(sum(track_A2B_prev{ii}>0));
    std_strain_23(ii,1) = std(E_total{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2))/2,[],'all','omitnan')/sqrt(sum(track_A2B_prev{ii}>0));
end

%u_total = inc2cum(u,dm,m,'linear');


%
step_num = 1:length(mean_strain_11);
figure
shadedErrorBar(step_num,mean_strain_11,std_strain_11,'b-x',1)
hold on
shadedErrorBar(step_num,mean_strain_22,std_strain_22,'g-*',1)
shadedErrorBar(step_num,mean_strain_33,std_strain_33,'r-+',1)
shadedErrorBar(step_num,mean_strain_12,std_strain_12,'m-o',1)
shadedErrorBar(step_num,mean_strain_13,std_strain_13,'y-s',1)
shadedErrorBar(step_num,mean_strain_23,std_strain_23,'k-^',1)
xlabel('step num')
ylabel('Strain')
axis([0,80,-.3,0.3])

title('Shear; std err shaded region')

