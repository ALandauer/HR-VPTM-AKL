
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use a scattered2grid and incremental finite deformation multiplicative cumulation
%

% set up vars
smoothness = MPTPara.smoothness;
grid_spacing = MPTPara.spacing; 


%define data range, pad edges to accond for bead near edges moving outside
%the range of the particles in the first image
coords_cur = resultDisp{1}.parCoordB;
xList = MPTPara.xRange(1)-2*MPTPara.edge_width*MPTPara.axesScale(1):grid_spacing(1):MPTPara.xRange(2)+2*MPTPara.edge_width*MPTPara.axesScale(1);
yList = MPTPara.yRange(1)-2*MPTPara.edge_width*MPTPara.axesScale(2):grid_spacing(2):MPTPara.yRange(2)+2*MPTPara.edge_width*MPTPara.axesScale(2);
zRange_grid(1) = min(coords_cur(:,3))-2*grid_spacing(3);
zRange_grid(2) = max(coords_cur(:,3))+2*grid_spacing(3);
zList = min(zRange_grid):grid_spacing(3):max(zRange_grid);
[yGrid,xGrid,zGrid] = meshgrid(yList,xList,zList);

%Set z range without padding
MPTPara.zRange(1) = min(coords_cur(:,3));
MPTPara.zRange(2) = max(coords_cur(:,3));

disp('%%%%% Interpolating tracking results on ref grid %%%%%'); fprintf('\n');
%interpolate scattered data onto eulerian grid
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
                   coords_cur_(:,3) < min(MPTPara.zRange)|...
                   coords_cur_(:,3) > max(MPTPara.zRange)];
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
                    coords_cur_(:,3) < min(MPTPara.zRange)|...
                    coords_cur_(:,3) > max(MPTPara.zRange)];
    coords_cur_strain = coords_cur_(keep_coords,:);
                                         
    
    F_cur_vec = resultDefGrad{ii}.F_A2B_refA;
    F_cur_ = zeros(length(coords_cur_),6);
    F_cur_(:,1) = F_cur_vec(1:9:end);
    F_cur_(:,2) = F_cur_vec(5:9:end);
    F_cur_(:,3) = F_cur_vec(9:9:end);
    F_cur_(:,4) = F_cur_vec(4:9:end);
    F_cur_(:,5) = F_cur_vec(7:9:end);
    F_cur_(:,6) = F_cur_vec(8:9:end);
    F_cur = F_cur_(keep_coords,:);
    
    [~,~,~,H_inc{ii}{1,1}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,1),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{2,2}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,2),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{3,3}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,3),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{1,2}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,4),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{1,3}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,5),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{2,3}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,6),grid_spacing,smoothness,xGrid,yGrid,zGrid);
    H_inc{ii}{2,1} = H_inc{ii}{1,2}; H_inc{ii}{3,1} = H_inc{ii}{1,3}; H_inc{ii}{3,2} = H_inc{ii}{2,3};
    
    [~,~,~,nan_mask_F{ii}{1}] = funScatter2Grid3D(coords_cur_strain(:,1),coords_cur_strain(:,2),coords_cur_strain(:,3),F_cur(:,1),grid_spacing,0,xGrid,yGrid,zGrid);
    nan_mask_F{ii}{1} = ones(size(nan_mask_F{ii}{1}))*nan_mask_F{ii}{1}./nan_mask_F{ii}{1};
    
end

%%
% compute the cumulative displacement field
disp('%%%%% Cumulating displacements %%%%%'); fprintf('\n');
msh{1} = xGrid;
msh{2} = yGrid;
msh{3} = zGrid;
[u_total,nan_mask] = inc2cum(u_inc,grid_spacing(1),msh,'linear');


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

% F_total_ = calculateFij(u_total,grid_spacing(1),[1,1,1],'optimal5');
% [E_total,e_total] = calculateEij(F_total_);

[E_total,e_total] = calculateEij(F_total);

%% save the results
results_file_names_Eul = fullfile('results',['results_3D_EulTotal_',file_names{1}(1:end-4),'.mat']);
if ~exist('results','dir') 
   mkdir('results')
end

resultDispInc = resultDisp;
resultDefGradInc = resultDefGrad;
save(results_file_names_Eul,'resultDispInc','resultDefGradInc','u_total','F_total','E_total','beadParam_all','MPTPara');

disp('%%%%% Cumulated Eulerian reuslts saved %%%%%'); fprintf('\n');

%%%%%% plotting options %%%%%%%
%slice views:

%% disp plots
clear H
cmr_colors = cmrMap(256);
frame_num = 10;

%adjust bd_wd and multipers to crop out edge effects
bd_wd = 4;
disp_vol_x = nan_mask{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*u_total{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
disp_vol_y = nan_mask{frame_num}{2}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*u_total{frame_num}{2}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
disp_vol_z = nan_mask{frame_num}{3}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*u_total{frame_num}{3}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));

yGrid_crop = yGrid(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
xGrid_crop = xGrid(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
zGrid_crop = zGrid(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));


figure
subplot(1,3,1)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,disp_vol_x,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-0.15,0.15])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('x-displacement')
colorbar

subplot(1,3,2)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,disp_vol_y,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-0.15,0.15])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('y-displacement')
colorbar

subplot(1,3,3)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,disp_vol_z,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-150,10])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('z-displacement')
colorbar

%% strain plots
clear H
cmr_colors = cmrMap(256);
frame_num = 10;

bd_wd = 4;
e_vol_xx = nan_mask_F{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*E_total{frame_num}{1,1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
e_vol_yy = nan_mask_F{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*E_total{frame_num}{2,2}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
e_vol_zz = nan_mask_F{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*E_total{frame_num}{3,3}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
e_vol_xy = nan_mask_F{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*E_total{frame_num}{1,2}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
e_vol_xz = nan_mask_F{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*E_total{frame_num}{1,3}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
e_vol_yz = nan_mask_F{frame_num}{1}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1)).*E_total{frame_num}{2,3}(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));

yGrid_crop = yGrid(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
xGrid_crop = xGrid(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));
zGrid_crop = zGrid(1*bd_wd:end-1*bd_wd,bd_wd:end-1*bd_wd,round(bd_wd*1):end-round(bd_wd*1));

xs = 700;
ys = 700;
zs = -300;

figure
subplot(2,3,1)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,e_vol_xx,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-0.15,0.15])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('E_{xx}')
colorbar

subplot(2,3,2)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,e_vol_yy,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-0.15,0.15])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('E_{yy}')
colorbar

subplot(2,3,3)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,e_vol_zz,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-0.15,0.30])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('E_{zz}')
colorbar

subplot(2,3,4)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,e_vol_xy,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-0.10,0.10])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('E_{xy}')
colorbar


subplot(2,3,5)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,e_vol_xz,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([0.25,0.35])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('E_{zx}')
colorbar

subplot(2,3,6)
H = slice(yGrid_crop,xGrid_crop,zGrid_crop,e_vol_yz,xs,ys,zs);
colormap(cmr_colors)
axis image
% caxis([-0.10,0.10])

H(1).EdgeColor = 'none';
H(2).EdgeColor = 'none';
H(3).EdgeColor = 'none';
xlabel('y (\mum)')
ylabel('x (\mum)')
zlabel('z (\mum)')
title('E_{zy}')
colorbar

