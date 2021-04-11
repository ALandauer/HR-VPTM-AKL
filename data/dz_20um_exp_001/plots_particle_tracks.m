parCoordTrajMat = cell2mat(parCoordTraj);
clear disp_A2BCum RMSD_* parCoordA parCoordB parCoordACum parCoordBCum parPosHist

[row1,col1] = find(isnan(parCoordTrajMat(1:length(file_names):end,1))==0);
trackParCum_ind = row1;
trackParCum_track_ratio = [];
%% z
figure
for ii = 1:20

for ImgSeqNum = 2:length(file_names)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(file_names):end,1))==0);
    trackParCum_ind = intersect(row2,trackParCum_ind);
    trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) / size(parCoord_prev{1},1);
    
    parCoordA = parCoordTrajMat(1:length(file_names):end,1:3);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(file_names):end,1:3);
    parCoordACum{ImgSeqNum-1} = parCoordA(trackParCum_ind,:);
    parCoordBCum{ImgSeqNum-1} = parCoordB(trackParCum_ind,:);
    disp_A2BCum{ImgSeqNum-1} = parCoordBCum{ImgSeqNum-1} - parCoordACum{ImgSeqNum-1};
    
    try
    parNum = find(trackParCum_ind == ii);
    parPosHist(ImgSeqNum-1) = parCoordBCum{ImgSeqNum-1}(parNum,3)'; 
    end
end


z_locs = [-180:-20:-620];

subplot(1,2,1)
hold on
plot(z_locs,-(parPosHist))
axis image
xlabel('z-height of particle plane')
ylabel('particle positions')
subplot(1,2,2)
hold on
plot(z_locs,abs(abs(parPosHist)-abs(z_locs)))
axis image
xlabel('z-height of particle plane')
ylabel('absolute position error')
end
legend('particle# 1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','19','20')




% % % parCoord_from_pisson_disk
% % % parCoord_from_localization
% % % 
% % % matches_A2B = f_track_nearest_neighbour3( parCoord_from_pisson_disk, parCoord_from_localization, Inf )
% % % [track_A2B, u_A2B] = funCompDisp3(parCoord_from_pisson_disk, parCoord_from_localization, matches_A2B, 0);





