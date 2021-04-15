clear
close all
load('C:\Users\akl1\Desktop\LFM\local\data\Tension_synthetic_volume\results_tension_complete_0414.mat')
init_config = load('C:\Users\akl1\Desktop\LFM\local\data\Tension_synthetic_volume\setup_stretch_25perc_syn.mat');

%%
disp('%%%%% Compute tracked cumulative displacements %%%%%'); fprintf('\n');
parCoordTrajMat = cell2mat(parCoordTraj);
clear disp_A2BCum* RMSD_* parCoordA parCoordB parCoordACum parCoordBCum parPosHist parCoordA_GT

[row1,col1] = find(isnan(parCoordTrajMat(1:length(file_names):end,1))==0);
trackParCum_ind = row1;
trackParCum_track_ratio = [];

deform_func = 1.05:.05:1.25; 
expt_strain22 = deform_func-1;
expt_strain11 = 1./sqrt(deform_func)-1;
expt_strain33 = 1./sqrt(deform_func)-1;

for ImgSeqNum = 2:length(file_names)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(file_names):end,1))==0);
    trackParCum_ind = intersect(row2,trackParCum_ind);
    trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) / size(parCoord_prev{1},1);
    
    parCoordA = parCoordTrajMat(1:length(file_names):end,1:3);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(file_names):end,1:3);
    
    parCoordACum{ImgSeqNum-1} = parCoordA(trackParCum_ind,:);
    parCoordBCum{ImgSeqNum-1} = parCoordB(trackParCum_ind,:);
    disp_A2BCum_exp{ImgSeqNum-1} = parCoordBCum{ImgSeqNum-1} - parCoordACum{ImgSeqNum-1};
    
    sizeI = [812,1288];%size(MPTPara.ImgRefMask);
    sizeV = [797,1275,161];%size(vol_in)
    parCoordA_GT{ImgSeqNum-1}(:,1) = parCoordACum{ImgSeqNum-1}(:,1) - MPTPara.axesScale(1)*sizeI(1)/2;
    parCoordA_GT{ImgSeqNum-1}(:,2) = parCoordACum{ImgSeqNum-1}(:,2);% - MPTPara.axesScale(2)*sizeI(2)/2;
    parCoordA_GT{ImgSeqNum-1}(:,3) = parCoordACum{ImgSeqNum-1}(:,3);% + MPTPara.axesScale(3)*sizeV(3)/2;

    dy = deform_func(ImgSeqNum-1);
    y0 = parCoordA_GT{ImgSeqNum-1};
    y1 = parCoordA_GT{ImgSeqNum-1};
    
    y1(:,1) = y1(:,1)/(sqrt(dy));
    y1(:,2) = y1(:,2)*dy;
    y1(:,3) = y1(:,3)/(sqrt(dy));
    parCoordBCum_GT{ImgSeqNum-1} = y1;
    
    disp_A2BCum_GT{ImgSeqNum-1} = y1 - y0;
    
    
    %     parNum = find(trackParCum_ind == 10);
    %     parPosHist(ImgSeqNum-1) = parCoordBCum{ImgSeqNum-1}(parNum,3)';
    
end

figure
quiver3(y0(:,1),y0(:,2),y0(:,3),disp_A2BCum_GT{end}(:,1),disp_A2BCum_GT{end}(:,2),disp_A2BCum_GT{end}(:,3))
hold on
quiver3(y0(:,1),y0(:,2),y0(:,3),disp_A2BCum_exp{end}(:,1),disp_A2BCum_exp{end}(:,2),disp_A2BCum_exp{end}(:,3))
legend('GT','Exp')

figure
plot3(parCoordACum{1}(:,2),parCoordACum{1}(:,1),parCoordACum{1}(:,3),'or')
hold on
plot3(init_config.LF_y0_ums{1}(:,1),init_config.LF_y0_ums{1}(:,2),init_config.LF_y0_ums{1}(:,3),'ob')
legend('GT','Exp')

% % mean_cum_disp = cellfun(@(x) mean(x,1),disp_A2BCum_exp,'UniformOutput',false);
% % mean_cum_disp = reshape(cell2mat(mean_cum_disp),3,[])';
% % 
% % std_cum_disp = cellfun(@(x) std(x,[],1),disp_A2BCum_exp,'UniformOutput',false);
% % std_cum_disp = reshape(cell2mat(std_cum_disp),3,[])';

%%
for ii = 1:length(disp_A2BCum_exp)
    disp_meas_x = disp_A2BCum_exp{ii}(:,2);
    disp_meas_y = disp_A2BCum_exp{ii}(:,1);
    disp_meas_z = disp_A2BCum_exp{ii}(:,3);
    
    N = length(disp_meas_z);
    
    disp_imps_x = disp_A2BCum_GT{ii}(:,2);
    disp_imps_y = disp_A2BCum_GT{ii}(:,1);
    disp_imps_z = disp_A2BCum_GT{ii}(:,3);
    
    RMSD_disp_y(ii,1) = sqrt(sum((disp_meas_y - disp_imps_y).^2)/N);
    RMSD_disp_x(ii,1) = sqrt(sum((disp_meas_x - disp_imps_x).^2)/N);
    RMSD_disp_z(ii,1) = sqrt(sum((disp_meas_z - disp_imps_z).^2)/N);
    
    mean_err_x(ii,1) = mean((disp_meas_x - disp_imps_x));
    mean_err_y(ii,1) = mean((disp_meas_y - disp_imps_y));
    mean_err_z(ii,1) = mean((disp_meas_z - disp_imps_z));
    
    std_err_x(ii,1) = std((disp_meas_x - disp_imps_x));
    std_err_y(ii,1) = std((disp_meas_y - disp_imps_y));
    std_err_z(ii,1) = std((disp_meas_z - disp_imps_z));
    
end

% disp(disp_meas_z - disp_imps_z)

x_lbl = 'Applied strain';

figure
shadedErrorBar(expt_strain22,mean_err_x,RMSD_disp_x,'g-x',1)
% shadedErrorBar(imps_disp,mean_cum_disp(:,2),std_cum_disp(:,2))
% xlabel('Noise level')
xlabel(x_lbl)
ylabel('Error level, um')
% axis image

hold on
shadedErrorBar(expt_strain22,mean_err_y,RMSD_disp_y,'b:*',1)
% shadedErrorBar(imps_disp,mean_cum_disp(:,1),std_cum_disp(:,1))
% axis image
shadedErrorBar(expt_strain22,mean_err_z,RMSD_disp_z,'r-.+',1)
% shadedErrorBar(imps_disp,mean_cum_disp(:,3),std_cum_disp(:,3))
title('RMSD shaded region')
% axis image

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

clear mean_strain_* std_strain_* RMSD_strain_* u_total

for ii = 1:length(disp_A2BCum_exp)
    
    disp_A2BCum_vec = disp_A2BCum_exp{ii};
    parCoordACum_vec = parCoordACum{ii};
    parCoordBCum_vec = parCoordBCum{ii};
%     parCoordBCum_vec = parCoordBCum_GT{ii};
    
    sxyz = min([round(0.5*MPTPara.f_o_s),20]);%.*MPTPara.axesScale; % Step size for griddata
    smoothness = 5e-2; % Smoothness for regularization; "smoothness=0" means no regularization
    
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
    
     disp_A2BCum_vec_GT = disp_A2BCum_GT{ii};
    parCoordACum_vec_GT = parCoordACum{ii};
    parCoordBCum_vec_GT = parCoordBCum_GT{ii};
    
    [~,~,~,u_Grid_refB_GT] = funScatter2Grid3D(parCoordACum_vec_GT(:,1),parCoordACum_vec_GT(:,2),parCoordACum_vec_GT(:,3),disp_A2BCum_vec_GT(:,1),sxyz,smoothness);
    [~,~,~,v_Grid_refB_GT] = funScatter2Grid3D(parCoordACum_vec_GT(:,1),parCoordACum_vec_GT(:,2),parCoordACum_vec_GT(:,3),disp_A2BCum_vec_GT(:,2),sxyz,smoothness);
    [~,~,~,w_Grid_refB_GT] = funScatter2Grid3D(parCoordACum_vec_GT(:,1),parCoordACum_vec_GT(:,2),parCoordACum_vec_GT(:,3),disp_A2BCum_vec_GT(:,3),sxyz,smoothness);
    
    % Build a displacement vector
    uvw_Grid_refB_Vector_GT = [u_Grid_refB_GT(:),v_Grid_refB_GT(:),w_Grid_refB_GT(:)]';
    uvw_Grid_refB_Vector_GT = uvw_Grid_refB_Vector_GT(:);
    % Calculate deformation gradient
    F_Grid_refB_Vector_GT = D_Grid*uvw_Grid_refB_Vector_GT; % {F} - I ={D}{U}
    
    %exp strain from def grad
    eij_exp{ii}{1,1} = reshape(F_Grid_refB_Vector(1:9:end),size(x_Grid_refB));
    eij_exp{ii}{2,2} = reshape(F_Grid_refB_Vector(5:9:end),size(x_Grid_refB));
    eij_exp{ii}{3,3} = reshape(F_Grid_refB_Vector(9:9:end),size(x_Grid_refB));
    eij_exp{ii}{1,2} = reshape(0.5*(F_Grid_refB_Vector(2:9:end)+F_Grid_refB_Vector(4:9:end)),size(x_Grid_refB));
    eij_exp{ii}{1,3} = reshape(0.5*(F_Grid_refB_Vector(3:9:end)+F_Grid_refB_Vector(7:9:end)),size(x_Grid_refB));
    eij_exp{ii}{2,3} = reshape(0.5*(F_Grid_refB_Vector(6:9:end)+F_Grid_refB_Vector(8:9:end)),size(x_Grid_refB));
    eij_exp{ii}{2,1} = eij_exp{ii}{1,2}; eij_exp{ii}{3,1} = eij_exp{ii}{1,3}; eij_exp{ii}{3,2} = eij_exp{ii}{2,3};
    
    %GT strain
    eij_GT{ii}{1,1} = reshape(F_Grid_refB_Vector_GT(1:9:end),size(x_Grid_refB));
    eij_GT{ii}{2,2} = reshape(F_Grid_refB_Vector_GT(5:9:end),size(x_Grid_refB));
    eij_GT{ii}{3,3} = reshape(F_Grid_refB_Vector(9:9:end),size(x_Grid_refB));
    eij_GT{ii}{1,2} = reshape(0.5*(F_Grid_refB_Vector_GT(2:9:end)+F_Grid_refB_Vector_GT(4:9:end)),size(x_Grid_refB));
    eij_GT{ii}{1,3} = reshape(0.5*(F_Grid_refB_Vector_GT(3:9:end)+F_Grid_refB_Vector_GT(7:9:end)),size(x_Grid_refB));
    eij_GT{ii}{2,3} = reshape(0.5*(F_Grid_refB_Vector_GT(6:9:end)+F_Grid_refB_Vector_GT(8:9:end)),size(x_Grid_refB));
    eij_GT{ii}{2,1} = eij_GT{ii}{1,2}; eij_GT{ii}{3,1} = eij_GT{ii}{1,3}; eij_GT{ii}{3,2} = eij_GT{ii}{2,3};
    
    mean_strain_err_11(ii,1) = mean(eij_exp{ii}{1,1}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,1}(11:end-10,11:end-10,3:end-2),'all','omitnan');
    mean_strain_err_22(ii,1) = mean(eij_exp{ii}{2,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,2}(11:end-10,11:end-10,3:end-2),'all','omitnan');
    mean_strain_err_33(ii,1) = mean(eij_exp{ii}{3,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{3,3}(11:end-10,11:end-10,3:end-2),'all','omitnan');
    mean_strain_err_12(ii,1) = mean(eij_exp{ii}{1,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,2}(11:end-10,11:end-10,3:end-2),'all','omitnan');
    mean_strain_err_13(ii,1) = mean(eij_exp{ii}{1,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,3}(11:end-10,11:end-10,3:end-2),'all','omitnan');
    mean_strain_err_23(ii,1) = mean(eij_exp{ii}{2,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,3}(11:end-10,11:end-10,3:end-2),'all','omitnan');
    
    mean_strain_13(ii,1) = mean(eij_exp{ii}{1,3}(:),'omitnan');
    mean_strain_23(ii,1) = mean(eij_exp{ii}{2,3}(:),'omitnan');
    
    std_strain_err_11(ii,1) = std(eij_exp{ii}{1,1}(:) - eij_GT{ii}{1,1}(:),'omitnan');
    std_strain_err_22(ii,1) = std(eij_exp{ii}{2,2}(:) - eij_GT{ii}{2,2}(:),'omitnan');
    std_strain_err_33(ii,1) = std(eij_exp{ii}{3,3}(:) - eij_GT{ii}{3,3}(:),'omitnan');
    std_strain_err_12(ii,1) = std(eij_exp{ii}{1,2}(:) - eij_GT{ii}{1,2}(:),'omitnan');
    std_strain_err_13(ii,1) = std(eij_exp{ii}{1,3}(:) - eij_GT{ii}{1,3}(:),'omitnan');
    std_strain_err_23(ii,1) = std(eij_exp{ii}{2,3}(:) - eij_GT{ii}{2,3}(:),'omitnan');
    
    N = numel(eij_exp{ii}{1,1});
    RMSD_strain_11(ii,1) = sqrt(sum((eij_exp{ii}{1,1}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,1}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
    RMSD_strain_22(ii,1) = sqrt(sum((eij_exp{ii}{2,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,2}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
    RMSD_strain_33(ii,1) = sqrt(sum((eij_exp{ii}{3,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{3,3}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
    RMSD_strain_12(ii,1) = sqrt(sum((eij_exp{ii}{1,2}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,2}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
    RMSD_strain_13(ii,1) = sqrt(sum((eij_exp{ii}{1,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{1,3}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
    RMSD_strain_23(ii,1) = sqrt(sum((eij_exp{ii}{2,3}(11:end-10,11:end-10,3:end-2) - eij_GT{ii}{2,3}(11:end-10,11:end-10,3:end-2)).^2,'all','omitnan')/N);
    
end
%%
figure
shadedErrorBar(expt_strain22,mean_strain_err_11,RMSD_strain_11,'g-x',1)
% shadedErrorBar(imps_disp,mean_cum_disp(:,2),std_cum_disp(:,2))
% xlabel('Noise level')
hold on
shadedErrorBar(expt_strain22,mean_strain_err_22,RMSD_strain_22,'b-*',1)
% shadedErrorBar(imps_disp,mean_cum_disp(:,1),std_cum_disp(:,1))
% axis image
shadedErrorBar(expt_strain22,mean_strain_err_33,RMSD_strain_33,'r-+',1)
% shadedErrorBar(imps_disp,mean_cum_disp(:,3),std_cum_disp(:,3))

shadedErrorBar(expt_strain22,mean_strain_err_12,RMSD_strain_12,'m-.o',1)
shadedErrorBar(expt_strain22,mean_strain_err_13,RMSD_strain_13,'y-.s',1)
shadedErrorBar(expt_strain22,mean_strain_err_23,RMSD_strain_23,'k-.^',1)
title('RMSD shaded region')
xlabel(x_lbl)
ylabel('Strain error')

title('Vol conserving tension; RMSD shaded region')