clear
close all
load('results_shear_exp_complete_35_to_100.mat')

%%
% disp('%%%%% Compute tracked cumulative displacements %%%%%'); fprintf('\n');
% parCoordTrajMat = cell2mat(parCoordTraj);
% clear disp_A2BCum* RMSD_disp* parCoordA parCoordB parCoordACum parCoordBCum parPosHist parCoordA_GT
% 
% [row1,col1] = find(isnan(parCoordTrajMat(1:length(file_names):end,1))==0);
% trackParCum_ind = row1;
% trackParCum_track_ratio = [];
% 
% for ImgSeqNum = 2:length(file_names)
%     [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(file_names):end,1))==0);
%     trackParCum_ind = intersect(row2,trackParCum_ind);
%     trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) / size(parCoord_prev{1},1);
%     
%     parCoordA = parCoordTrajMat(1:length(file_names):end,1:3);
%     parCoordB = parCoordTrajMat(ImgSeqNum:length(file_names):end,1:3);
%     
%     parCoordACum{ImgSeqNum-1} = parCoordA(trackParCum_ind,:);
%     parCoordBCum{ImgSeqNum-1} = parCoordB(trackParCum_ind,:);
%     disp_A2BCum_exp{ImgSeqNum-1} = parCoordBCum{ImgSeqNum-1} - parCoordACum{ImgSeqNum-1};
%     
% end

%%
clear u_inc F_ij* mean_strain_* std_strain_* H_inc*
u_inc = cell(1,3);
dm = 25;
sxyz = [dm,dm,dm];
smoothness = 0;

coords_cur = resultDisp{1}.parCoordA;
xList = 20:sxyz(1):1300+sxyz(1);
yList = 20:sxyz(2):1050+sxyz(2);
zList = -520:sxyz(3):100+sxyz(3);
[yGrid,xGrid,zGrid] = meshgrid(yList,xList,zList);

for ii = 1:5%length(uvw_B2A_prev)
    disp(ii)
    
    coords_cur = resultDisp{ii}.parCoordB;
    disps_cur = resultDisp{ii}.disp_A2B_parCoordB;
    [~,~,~,u_inc{ii}{1}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),disps_cur(:,1),sxyz,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,u_inc{ii}{2}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),disps_cur(:,2),sxyz,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,u_inc{ii}{3}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),disps_cur(:,3),sxyz,smoothness,xGrid,yGrid,zGrid);
    
    F_cur = resultDefGrad{ii}.F_A2B_refA;
    [~,~,~,H_inc{ii}{1,1}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),F_cur(1:9:end),sxyz,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{2,2}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),F_cur(2:9:end),sxyz,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{3,3}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),F_cur(3:9:end),sxyz,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{1,2}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),F_cur(4:9:end),sxyz,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{1,3}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),F_cur(6:9:end),sxyz,smoothness,xGrid,yGrid,zGrid);
    [~,~,~,H_inc{ii}{2,3}] = funScatter2Grid3D(coords_cur(:,1),coords_cur(:,2),coords_cur(:,3),F_cur(8:9:end),sxyz,smoothness,xGrid,yGrid,zGrid);
    H_inc{ii}{2,1} = H_inc{ii}{1,2}; H_inc{ii}{3,1} = H_inc{ii}{1,3}; H_inc{ii}{3,2} = H_inc{ii}{2,3};
end

%%
m{1} = xGrid;
m{2} = yGrid;
m{3} = zGrid;
u_total = inc2cum(u_inc,dm,m,'linear');


%%
F_total = cell(length(H_inc)+1,1);
F_total{1}{1,1} = ones(size(H_inc{1}{1,1}));
F_total{1}{2,2} = ones(size(H_inc{1}{2,2}));
F_total{1}{3,3} = ones(size(H_inc{1}{3,3}));

F_total{1}{1,2} = zeros(size(H_inc{1}{1,2}));
F_total{1}{2,1} = zeros(size(H_inc{1}{2,1}));
F_total{1}{1,3} = zeros(size(H_inc{1}{1,3}));
F_total{1}{3,1} = zeros(size(H_inc{1}{3,1}));
F_total{1}{2,3} = zeros(size(H_inc{1}{2,3}));
F_total{1}{3,2} = zeros(size(H_inc{1}{3,2}));

for ii = 1:length(H_inc)
    disp(ii)
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
            F_total{ii+1}{1,1}(i,j,k) = Fhat_inc(1,1,loc);
            F_total{ii+1}{2,2}(i,j,k) = Fhat_inc(2,2,loc);
            F_total{ii+1}{3,3}(i,j,k) = Fhat_inc(3,3,loc);
            F_total{ii+1}{1,2}(i,j,k) = Fhat_inc(1,2,loc);
            F_total{ii+1}{2,1}(i,j,k) = Fhat_inc(2,1,loc);
            F_total{ii+1}{1,3}(i,j,k) = Fhat_inc(1,3,loc);
            F_total{ii+1}{3,1}(i,j,k) = Fhat_inc(3,1,loc);
            F_total{ii+1}{2,3}(i,j,k) = Fhat_inc(2,3,loc);
            F_total{ii+1}{3,2}(i,j,k) = Fhat_inc(2,3,loc);
            
        else
            
            Fcur(1,1,loc) = F_total{ii}{1,1}(loc);
            Fcur(2,2,loc) = F_total{ii}{2,2}(loc);
            Fcur(3,3,loc) = F_total{ii}{3,3}(loc);
            
            Fcur(1,2,loc) = F_total{ii}{1,2}(loc);
            Fcur(2,1,loc) = F_total{ii}{2,1}(loc);
            Fcur(1,3,loc) = F_total{ii}{1,3}(loc);
            Fcur(3,1,loc) = F_total{ii}{3,1}(loc);
            Fcur(2,3,loc) = F_total{ii}{2,3}(loc);
            Fcur(3,2,loc) = F_total{ii}{3,2}(loc);
            
            
            %F_total = F0*F1*F2...Fn
            F_tot_ = Fcur(:,:,loc)*Fhat_inc(:,:,loc);
            
            F_total{ii+1}{1,1}(i,j,k) = F_tot_(1,1);
            F_total{ii+1}{2,2}(i,j,k) = F_tot_(2,2);
            F_total{ii+1}{3,3}(i,j,k) = F_tot_(3,3);
            
            F_total{ii+1}{1,2}(i,j,k) = F_tot_(1,2);
            F_total{ii+1}{2,1}(i,j,k) = F_tot_(2,1);
            F_total{ii+1}{1,3}(i,j,k) = F_tot_(1,3);
            F_total{ii+1}{3,1}(i,j,k) = F_tot_(3,1);
            F_total{ii+1}{2,3}(i,j,k) = F_tot_(2,3);
            F_total{ii+1}{3,2}(i,j,k) = F_tot_(3,2);
            
        end
    end
    
end
%%
for ii = 1:length(F_total)
    bd_wd = 5;
    
    mean_strain_11(ii,1) = mean(F_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan')-1;
    mean_strain_22(ii,1) = mean(F_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan')-1;
    mean_strain_33(ii,1) = mean(F_total{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan')-1;
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
shadedErrorBar(step_num,mean_strain_13,std_strain_13,'m-.s',1)
shadedErrorBar(step_num,mean_strain_23,std_strain_23,'k-.^',1)
xlabel('step num')
ylabel('Strain')

title('Shear; Stdev shaded region')
%%
% % % bd_wd = 5;
% % %
% % % mean_strain_11(ii,1) = mean(eij_exp{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % % mean_strain_22(ii,1) = mean(eij_exp{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % % mean_strain_33(ii,1) = mean(eij_exp{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % % mean_strain_12(ii,1) = mean(eij_exp{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % % mean_strain_13(ii,1) = mean(eij_exp{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % % mean_strain_23(ii,1) = mean(eij_exp{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % %
% % % %     mean_strain_13(ii,1) = mean(eij_exp{ii}{1,3}(:),'omitnan');
% % % %     mean_strain_23(ii,1) = mean(eij_exp{ii}{2,3}(:),'omitnan');
% % %
% % % std_strain_11(ii,1) = std(eij_exp{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
% % % std_strain_22(ii,1) = std(eij_exp{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
% % % std_strain_33(ii,1) = std(eij_exp{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
% % % std_strain_12(ii,1) = std(eij_exp{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
% % % std_strain_13(ii,1) = std(eij_exp{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');
% % % std_strain_23(ii,1) = std(eij_exp{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),[],'all','omitnan');

% for ii = 1:length(disp_A2BCum_exp)
%     disp_meas_x = disp_A2BCum_exp{ii}(:,2);
%     disp_meas_y = disp_A2BCum_exp{ii}(:,1);
%     disp_meas_z = disp_A2BCum_exp{ii}(:,3);
%
%     N = length(disp_meas_z);
%
%     disp_imps_x = disp_A2BCum_GT{ii}(:,2);
%     disp_imps_y = disp_A2BCum_GT{ii}(:,1);
%     disp_imps_z = disp_A2BCum_GT{ii}(:,3);
%
%     RMSD_disp_x(ii,1) = sqrt(sum((disp_meas_x - disp_imps_x).^2)/N);
%     RMSD_disp_y(ii,1) = sqrt(sum((disp_meas_y - disp_imps_y).^2)/N);
%     RMSD_disp_z(ii,1) = sqrt(sum((disp_meas_z - disp_imps_z).^2)/N);
%
%     mean_disp_err_x(ii,1) = mean((disp_meas_x - disp_imps_x));
%     mean_disp_err_y(ii,1) = mean((disp_meas_y - disp_imps_y));
%     mean_disp_err_z(ii,1) = mean((disp_meas_z - disp_imps_z));
%
%     std_disp_err_x(ii,1) = std((disp_meas_x - disp_imps_x));
%     std_disp_err_y(ii,1) = std((disp_meas_y - disp_imps_y));
%     std_disp_err_z(ii,1) = std((disp_meas_z - disp_imps_z));
%
% end

% disp(disp_meas_z - disp_imps_z)

% x_lbl = 'Applied strain';
%
% figure
% shadedErrorBar(-deform_func,mean_disp_err_x,RMSD_disp_x,'g-x',1)
% % shadedErrorBar(imps_disp,mean_cum_disp(:,2),std_cum_disp(:,2))
% % xlabel('Noise level')
% xlabel(x_lbl)
% ylabel('Error level, um')
% % axis image
%
% hold on
% shadedErrorBar(-deform_func,mean_disp_err_y,RMSD_disp_y,'b:*',1)
% % shadedErrorBar(imps_disp,mean_cum_disp(:,1),std_cum_disp(:,1))
% % axis image
% shadedErrorBar(-deform_func,mean_disp_err_z,RMSD_disp_z,'r-.+',1)
% % shadedErrorBar(imps_disp,mean_cum_disp(:,3),std_cum_disp(:,3))
% title('RMSD shaded region')
% % axis image

%%

% % y1 = rotation(parCoordBCum{ii}, thetas(ii));
% % for ii = 2:length(thetas)
% %
% %     y1 = rotation(parCoordBCum{ii}, thetas(ii));
% %
% %     disp_A2BCumGT =
% % %     rot_coords(ii).rad = thetas(ii);
% % %     rot_coords(ii).deg = rad2deg(thetas(ii));
% %
% %     LF_y0_ums{ii}= y1;
% %     %rot_coords(ii).coords = y1;
% % end

%%

% % % clear mean_strain_* std_strain_* RMSD_strain_* u_total
% % % 
% % % for ii = 1:length(disp_A2BCum_exp)
% % %     disp(ii)
% % %     disp_A2BCum_vec = disp_A2BCum_exp{ii};
% % %     parCoordACum_vec = parCoordACum{ii};
% % %     parCoordBCum_vec = parCoordBCum{ii};
% % %     %     parCoordBCum_vec = parCoordBCum_GT{ii};
% % %     
% % %     sxyz = min([round(0.5*MPTPara.f_o_s),1]).*[20,20,10];%MPTPara.axesScale; % Step size for griddata
% % %     smoothness = 0.05; % Smoothness for regularization; "smoothness=0" means no regularization
% % %     
% % %     [x_Grid_refB,y_Grid_refB,z_Grid_refB,u_Grid_refB] = ...
% % %         funScatter2Grid3D(parCoordACum_vec(:,1),parCoordACum_vec(:,2),parCoordACum_vec(:,3),disp_A2BCum_vec(:,1),sxyz,smoothness);
% % %     [~,~,~,v_Grid_refB] = funScatter2Grid3D(parCoordACum_vec(:,1),parCoordACum_vec(:,2),parCoordACum_vec(:,3),disp_A2BCum_vec(:,2),sxyz,smoothness);
% % %     [~,~,~,w_Grid_refB] = funScatter2Grid3D(parCoordACum_vec(:,1),parCoordACum_vec(:,2),parCoordACum_vec(:,3),disp_A2BCum_vec(:,3),sxyz,smoothness);
% % %     
% % %     % Build a displacement vector
% % %     uvw_Grid_refB_Vector = [u_Grid_refB(:),v_Grid_refB(:),w_Grid_refB(:)]';
% % %     uvw_Grid_refB_Vector = uvw_Grid_refB_Vector(:);
% % %     % Calculate deformation gradient
% % %     D_Grid = funDerivativeOp3(size(x_Grid_refB,1),size(x_Grid_refB,2),size(x_Grid_refB,3), sxyz ); % Central finite difference operator
% % %     F_Grid_refB_Vector = D_Grid*uvw_Grid_refB_Vector; % {F} - I ={D}{U}
% % %     
% % %     %exp strain from def grad
% % %     eij_exp{ii}{1,1} = reshape(F_Grid_refB_Vector(1:9:end),size(x_Grid_refB));
% % %     eij_exp{ii}{2,2} = reshape(F_Grid_refB_Vector(5:9:end),size(x_Grid_refB));
% % %     eij_exp{ii}{3,3} = reshape(F_Grid_refB_Vector(9:9:end),size(x_Grid_refB));
% % %     eij_exp{ii}{1,2} = 2*reshape(0.5*(F_Grid_refB_Vector(2:9:end)+F_Grid_refB_Vector(4:9:end)),size(x_Grid_refB));
% % %     eij_exp{ii}{1,3} = 2*reshape(0.5*(F_Grid_refB_Vector(3:9:end)+F_Grid_refB_Vector(7:9:end)),size(x_Grid_refB));
% % %     eij_exp{ii}{2,3} = 2*reshape(0.5*(F_Grid_refB_Vector(6:9:end)+F_Grid_refB_Vector(8:9:end)),size(x_Grid_refB));
% % %     eij_exp{ii}{2,1} = eij_exp{ii}{1,2}; eij_exp{ii}{3,1} = eij_exp{ii}{1,3}; eij_exp{ii}{3,2} = eij_exp{ii}{2,3};
% % %     
% % %     bd_wd = 5;
% % %     
% % %     mean_strain_11(ii,1) = mean(eij_exp{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % %     mean_strain_22(ii,1) = mean(eij_exp{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % %     mean_strain_33(ii,1) = mean(eij_exp{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % %     mean_strain_12(ii,1) = mean(eij_exp{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % %     mean_strain_13(ii,1) = mean(eij_exp{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % %     mean_strain_23(ii,1) = mean(eij_exp{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/2):end-round(bd_wd/2)),'all','omitnan');
% % %     
% % %     %     mean_strain_13(ii,1) = mean(eij_exp{ii}{1,3}(:),'omitnan');
% % %     %     mean_strain_23(ii,1) = mean(eij_exp{ii}{2,3}(:),'omitnan');
% % %     
% % %     std_strain_11(ii,1) = std(eij_exp{ii}{1,1}(:),'omitnan');
% % %     std_strain_22(ii,1) = std(eij_exp{ii}{2,2}(:),'omitnan');
% % %     std_strain_33(ii,1) = std(eij_exp{ii}{3,3}(:),'omitnan');
% % %     std_strain_12(ii,1) = std(eij_exp{ii}{1,2}(:),'omitnan');
% % %     std_strain_13(ii,1) = std(eij_exp{ii}{1,3}(:),'omitnan');
% % %     std_strain_23(ii,1) = std(eij_exp{ii}{2,3}(:),'omitnan');
% % %     
% % %     %     mean_strain_err_11(ii,1) = mean(eij_exp{ii}{1,1}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,1}(11:end-10,11:end-10,3:end-2),'all','omitnan');
% % %     %     mean_strain_err_22(ii,1) = mean(eij_exp{ii}{2,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,2}(11:end-10,11:end-10,3:end-2),'all','omitnan');
% % %     %     mean_strain_err_33(ii,1) = mean(eij_exp{ii}{3,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{3,3}(11:end-10,11:end-10,3:end-2),'all','omitnan');
% % %     %     mean_strain_err_12(ii,1) = mean(eij_exp{ii}{1,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,2}(11:end-10,11:end-10,3:end-2),'all','omitnan');
% % %     %     mean_strain_err_13(ii,1) = mean(eij_exp{ii}{1,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,3}(11:end-10,11:end-10,3:end-2),'all','omitnan');
% % %     %     mean_strain_err_23(ii,1) = mean(eij_exp{ii}{2,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,3}(11:end-10,11:end-10,3:end-2),'all','omitnan');
% % %     %
% % %     %     mean_strain_13(ii,1) = mean(eij_exp{ii}{1,3}(:),'omitnan');
% % %     %     mean_strain_23(ii,1) = mean(eij_exp{ii}{2,3}(:),'omitnan');
% % %     %
% % %     %     std_strain_err_11(ii,1) = std(eij_exp{ii}{1,1}(:) - eij_GT{ii}{1,1}(:),'omitnan');
% % %     %     std_strain_err_22(ii,1) = std(eij_exp{ii}{2,2}(:) - eij_GT{ii}{2,2}(:),'omitnan');
% % %     %     std_strain_err_33(ii,1) = std(eij_exp{ii}{3,3}(:) - eij_GT{ii}{3,3}(:),'omitnan');
% % %     %     std_strain_err_12(ii,1) = std(eij_exp{ii}{1,2}(:) - eij_GT{ii}{1,2}(:),'omitnan');
% % %     %     std_strain_err_13(ii,1) = std(eij_exp{ii}{1,3}(:) - eij_GT{ii}{1,3}(:),'omitnan');
% % %     %     std_strain_err_23(ii,1) = std(eij_exp{ii}{2,3}(:) - eij_GT{ii}{2,3}(:),'omitnan');
% % %     
% % %     %     N = numel(eij_exp{ii}{1,1});
% % %     %     RMSD_strain_11(ii,1) = sqrt(sum((eij_exp{ii}{1,1}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,1}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
% % %     %     RMSD_strain_22(ii,1) = sqrt(sum((eij_exp{ii}{2,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,2}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
% % %     %     RMSD_strain_33(ii,1) = sqrt(sum((eij_exp{ii}{3,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{3,3}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
% % %     %     RMSD_strain_12(ii,1) = sqrt(sum((eij_exp{ii}{1,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,2}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
% % %     %     RMSD_strain_13(ii,1) = sqrt(sum((eij_exp{ii}{1,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,3}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
% % %     %     RMSD_strain_23(ii,1) = sqrt(sum((eij_exp{ii}{2,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,3}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
% % %     
% % % end
% % % %%
% % % step_num = 1:length(mean_strain_11);
% % % figure
% % % shadedErrorBar(step_num,mean_strain_11,std_strain_11,'g-x',1)
% % % % shadedErrorBar(imps_disp,mean_cum_disp(:,2),std_cum_disp(:,2))
% % % % xlabel('Noise level')
% % % hold on
% % % shadedErrorBar(step_num,mean_strain_22,std_strain_22,'b-*',1)
% % % % shadedErrorBar(imps_disp,mean_cum_disp(:,1),std_cum_disp(:,1))
% % % % axis image
% % % shadedErrorBar(step_num,mean_strain_33,std_strain_33,'r-+',1)
% % % % shadedErrorBar(imps_disp,mean_cum_disp(:,3),std_cum_disp(:,3))
% % % 
% % % shadedErrorBar(step_num,mean_strain_12,std_strain_12,'y-.o',1)
% % % shadedErrorBar(step_num,mean_strain_13,std_strain_13,'m-.s',1)
% % % shadedErrorBar(step_num,mean_strain_23,std_strain_23,'k-.^',1)
% % % xlabel('step num')
% % % ylabel('Strain')
% % % 
% % % title('Shear; Stdev shaded region')