%%% Compute cumulative displacement at each step %%%%%
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
    
    parNum = find(trackParCum_ind == 10);
    parPosHist(ImgSeqNum-1) = parCoordBCum{ImgSeqNum-1}(parNum,3)'; 
    
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
    
    sterr_y(ii,1) = std_cum_disp(ii,1)/sqrt(N);
    sterr_x(ii,1) = std_cum_disp(ii,2)/sqrt(N);
    sterr_z(ii,1) = std_cum_disp(ii,3)/sqrt(N);
end

x_lbl = 'Experimental imposed disp z';
% imps_disp = 20*[1:length(mean_cum_disp)]';
imps_disp = -20*[1:length(mean_cum_disp)]';
% imps_disp = [0.022,0.025,0.028,0.033,0.040,0.050,0.066,0.100,0.200];
figure
shadedErrorBar(imps_disp,mean_cum_disp(:,2),sterr_x,'b-x',1)
hold on
shadedErrorBar(imps_disp,mean_cum_disp(:,1),sterr_y,'g-*',1)
shadedErrorBar(imps_disp,mean_cum_disp(:,3),sterr_z,'r-+',1)
% shadedErrorBar(imps_disp,mean_cum_disp(:,3),std_cum_disp(:,3))
hold on
plot(imps_disp,imps_disp,'b--')
xlabel(x_lbl)
ylabel('Measured displacement in z, um')
title('St Err shaded region')
% axis image
