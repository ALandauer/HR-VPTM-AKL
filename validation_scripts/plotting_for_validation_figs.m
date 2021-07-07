%%
%plot mean strain comps (lagrangian tracking)
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

%%
%plot mean strain comps (eulerian)

%%
%plot mean strain comps
bd_wd = 4;
clear mean_strain_* std_strain_*


for ii = 1:length(E_total)-1

    N = sum(track_A2B_prev{ii}>0);
    
    mean_strain_11(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),'all','omitnan');
    mean_strain_22(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),'all','omitnan');
    mean_strain_33(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),'all','omitnan');
    mean_strain_12(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),'all','omitnan');
    mean_strain_13(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),'all','omitnan');
    mean_strain_23(ii,1) = mean(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),'all','omitnan');
    
    ste_strain_11(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{1,1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),[],'all','omitnan')/sqrt(N);
    ste_strain_22(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{2,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),[],'all','omitnan')/sqrt(N);
    ste_strain_33(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{3,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),[],'all','omitnan')/sqrt(N);
    ste_strain_12(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{1,2}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),[],'all','omitnan')/sqrt(N);
    ste_strain_13(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{1,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),[],'all','omitnan')/sqrt(N);
    ste_strain_23(ii,1) = std(nan_mask_F{ii}{1}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1))...
        .*E_total{ii}{2,3}(bd_wd:end-bd_wd,bd_wd:end-bd_wd,round(bd_wd/1):end-round(bd_wd/1)),[],'all','omitnan')/sqrt(N);
end

%u_total = inc2cum(u,dm,m,'linear');


%
step_num = 0.002*(0:length(mean_strain_11)-1);
figure
shadedErrorBar(step_num,mean_strain_11,ste_strain_11,'b-x',1)
hold on
shadedErrorBar(step_num,mean_strain_22,ste_strain_22,'g-*',1)
shadedErrorBar(step_num,mean_strain_33,ste_strain_33,'r-+',1)
shadedErrorBar(step_num,mean_strain_12,ste_strain_12,'m-o',1)
shadedErrorBar(step_num,mean_strain_13,ste_strain_13,'y-s',1)
shadedErrorBar(step_num,mean_strain_23,ste_strain_23,'k-^',1)
xlabel('Time, s')
ylabel('Lagrange strain')
% axis([0.002*0,0.002*85,-.1,0.25])

title('Shear; std err shaded region')


