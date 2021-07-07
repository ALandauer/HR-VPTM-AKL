%% Synthetic Image Generation
clear all; close;

%% Define volume and create centroids
sizeI = [1280,800,60]; %need to figure out depth of field issue and um2px conversion? FOVx, FOVy
spacing = 8;
nPts = 60;
y0 = poissonDisc(sizeI,spacing,nPts,0);

sigma = [1,1,1];

% Generate reference image
L0 = seedBeadsN(sigma,y0,sizeI);
%Figure volume
%imshow3D(L0,[]); 
imagesc3D(L0)
%%
%% Stretch with lambda in e1 direction

% for jj = 1:length(stretches)
%     y2=stretch(y0,stretches(jj));
%     stretch_coords(jj).stretch = stretches(jj);
%     stretch_coords(jj).coords = y2;
% end
% save('stretch_coords.mat','stretch_coords')

stretch = 1:.05:1.25; 
dy = stretch;

LF_y0(:,1) = y0(:,1)- sizeI(1)%/2;
LF_y0(:,2) = y0(:,2)- sizeI(2)/2;
LF_y0(:,3) = y0(:,3)- sizeI(3)/2;

LF_y0_um(:,1) = LF_y0(:,1)*1.64;
LF_y0_um(:,2) = LF_y0(:,2)*1.64;
LF_y0_um(:,3) = (LF_y0(:,3)*5)+310;


LF_ = LF_y0_um;%LF_y0_um;
for ii = 1:length(dy)
    
    LF(:,1) = LF_(:,1)*dy(ii);
    LF(:,2) = LF_(:,2)/(sqrt(dy(ii)));
    LF(:,3) = LF_(:,3)/(sqrt(dy(ii)));
    %LF_pixel{ii} = LF;
    
    LF(:,1) = -LF(:,1)- ((sizeI(1)/2) *1.64);
    LF(:,2) = LF(:,2)%- ((sizeI(2)/2) *1.64);
    LF(:,3) = LF(:,3)%- ((sizeI(3)/2) *5)
    
    LF_y0_ums{ii}= LF;
    
%     
%     % Convert pixels to LF-image space
%     LF_y0(:,1) = LF(:,1)- sizeI(1)/2
%     LF_y0(:,2) = LF(:,2)- sizeI(2)/2
%     LF_y0(:,3) = LF(:,3)- sizeI(3)/2 
%     
%     LF_y0_um(:,1)= LF_y0(:,1)*1.64;
%     LF_y0_um(:,2)= LF_y0(:,2)*1.64;
%     LF_y0_um(:,3)= LF_y0(:,3)*5 + 320;
%     LF_y0_all{ii} = LF_y0;
%     LF_y0_ums{ii} = LF_y0_um;
end



save('stretch_25perc_syn.mat','y0','sizeI','LF_y0','LF_y0_ums','nPts','stretch')


%% SIMPLE SHEAR CASE
per_strain = 0:.025:.25; %microns
dy = per_strain;
LF_ = y0;%LF_y0_um;

% Convert pixels to LF-image space
LF_y0(:,1) = y0(:,1)%- sizeI(1)/2;
LF_y0(:,2) = y0(:,2)%- sizeI(2)/2;
LF_y0(:,3) = y0(:,3)%- sizeI(3)/2;

LF_y0_um(:,1) = LF_y0(:,1)*1.64;
LF_y0_um(:,2) = LF_y0(:,2)*1.64;
LF_y0_um(:,3) = (LF_y0(:,3)*5)+310;

for ii = 1:length(dy)
    simp_shear(ii).perc = dy(ii);
    
    LF(:,2) = LF_y0_um(:,2)+(LF_y0_um(:,3)*dy(ii));
    LF(:,1) = LF_y0_um(:,1);
    LF(:,3) = LF_y0_um(:,3);
    
    LF(:,1) = LF(:,1)- ((sizeI(1)/2) *1.64);
    LF(:,2) = LF(:,2)- ((sizeI(2)/2) *1.64);
    LF(:,3) = LF(:,3)%- ((sizeI(3)/2) *5)%+310;


    %LF_pixel{ii} = LF;
    %LF_y0_all{ii} = LF_y0;
    
%     LF(:,1) = (LF(:,1)- (sizeI(1)/2)) *1.64;
%     LF(:,2) = (LF(:,2)- (sizeI(2)/2)) *1.64;
%     LF(:,3) = (LF(:,3)- (sizeI(3)/2)) *5 +310;
    
    LF_y0_ums{ii} = LF;
end
% LF_y0_ums(:,1)= LF(:,1)*1.64;
% LF_y0_ums(:,2)= LF(:,2)*1.64;
% LF_y0_ums(:,3)= (LF(:,3)*5)+310;

per_strain =dy ;
LF_y0_ums;

save('simple_shear_25perc_syn.mat','y0','sizeI','LF_y0','LF_y0_ums','nPts','per_strain')


%% RIGID BODY DISP OF SINGLE LAYER

% Convert pixels to LF-image space
LF_y0(:,1) = y0(:,1)- sizeI(1)/2
LF_y0(:,2) = y0(:,2)- sizeI(2)/2
LF_y0(:,3) = y0(:,3)- sizeI(3)/2 

LF_y0_um(:,1)= LF_y0(:,1)*1.64;
LF_y0_um(:,2)= LF_y0(:,2)*1.64;
LF_y0_um(:,3)= LF_y0(:,3)*5;
LF_y0_um(:,3) = LF_y0_um(:,3);


dy = 0:5:55; %microns
LF_ = LF_y0_um;
for ii = 1:length(dy)
    RBD(ii).dy = dy(ii);
    
    LF(:,1) = LF_y0_um(:,1)
    LF(:,2)=LF_y0_um(:,2)+dy(ii)
    LF(:,3)=LF_y0_um(:,3)+310;
    LF_y0_ums{ii} = LF;
end

%dz = dy;

% noise = [.02 .022 .025 .028 (1/30) (1/25) (1/20) (1/15) (1/10) (1/5)];
%LF_y0_ums = {LF_y0_um};
% save('syn_vol_noise.mat','y0','sizeI','LF_y0','LF_y0_um','nPts','noise2')

save('dy_5um_syn_vol.mat','y0','sizeI','LF_y0','LF_y0_um','LF_y0_ums','nPts','dy')
%% Rigid rotation through angle theta about e3 axis

LF_y0(:,1) = y0(:,1)- sizeI(1)/2
LF_y0(:,2) = y0(:,2)- sizeI(2)/2
LF_y0(:,3) = y0(:,3)- sizeI(3)/2 

LF_y0_um(:,1)= LF_y0(:,1)*1.64;
LF_y0_um(:,2)= LF_y0(:,2)*1.64;
LF_y0_um(:,3)= LF_y0(:,3)*5;
LF_y0_um(:,3) = LF_y0_um(:,3)+310;

clear thetas
thetas = deg2rad([0:5:45]); 
for ii = 1:length(thetas)
    
    y1=rotation(LF_y0_um, thetas(ii));
%     rot_coords(ii).rad = thetas(ii);
%     rot_coords(ii).deg = rad2deg(thetas(ii));
     
    LF_y0_ums{ii}= y1; 
    %rot_coords(ii).coords = y1;
end


save('rot_coords.mat','y0','sizeI','LF_y0','LF_y0_um','LF_y0_ums','nPts','thetas')
%%
% % kNN method on rotation validation
% for k = 1:length(rot_coords)-1
%     Idx = knnsearch(rot_coords(k+1).coords, rot_coords(k).coords);
%     u{k} = (rot_coords(k).coords - rot_coords(k+1).coords(Idx,:)); 
% end
% 
% % Analytical displacement
% for mm = 1:length(thetas)-1
%     u_analyt{mm}(:,1) = rot_coords(mm+1).coords(:,1)- rot_coords(mm).coords(:,1);
%     u_analyt{mm}(:,2) = rot_coords(mm+1).coords(:,2)- rot_coords(mm).coords(:,2);
%     u_analyt{mm}(:,3) =  rot_coords(mm+1).coords(:,3)- rot_coords(mm).coords(:,3);
% end
% 
% %RMS Error
% for kk = 1:length(u)
%     temp1 = u{kk} - u_analyt{kk};
%     for zz = 1:length(temp1)
%         sums = 0;
%         temp2 = (norm(temp1(zz,:)))^2;
%         sums = sums + temp2;
%     end
%     err{kk} = sqrt(sums)/length(y0);
% end
% 
% figure; 
% plot(1:length(u_analyt),cell2mat(err),'o')
% title('RMS Error wrt Rotation using kNN Search')
% xlabel('Applied rotation (deg)')
% ylabel('RMS')
% ylim([0 1])
% 
% %Quiver plots
% figure;
% for jj = 1:length(u)
%     quiver(rot_coords(jj).coords(:,1), rot_coords(jj).coords(:,2),u{jj}(:,1),u{jj}(:,2),0)
%     axis image
%     hold on
% end
% %% Generate Light-field images for rotation
% lfm_rot = LFM_RayTrace(1280, 800,38,1,.1,13,-4,605e-6,500000, y0);
% % imwrite(image8bit,'test_noise222.tif');
% % figure; imshow(lfm_rot,[])
% %% Stretch with lambda in e1 direction
% stretches = [0: .01:.1];
% for jj = 1:length(stretches)
%     y2=stretch(y0,stretches(jj));
%     stretch_coords(jj).stretch = stretches(jj);
%     stretch_coords(jj).coords = y2;
% end
% save('stretch_coords.mat','stretch_coords')
% %%
% % kNN method on rotation validation
% for k2 = 1:length(stretch_coords)-1
%     Idx = knnsearch(stretch_coords(k2+1).coords, stretch_coords(k2).coords);
%     u2{k2} = (stretch_coords(k2).coords - stretch_coords(k2+1).coords(Idx,:)); 
% end
% 
% % Analytical displacement
% for mm = 1:length(stretches)-1
%     u_analyt2{mm}(:,1) = stretch_coords(mm+1).coords(:,1)- stretch_coords(mm).coords(:,1);
%     u_analyt2{mm}(:,2) = stretch_coords(mm+1).coords(:,2)- stretch_coords(mm).coords(:,2);
%     u_analyt2{mm}(:,3) =  stretch_coords(mm+1).coords(:,3)- stretch_coords(mm).coords(:,3);
% end
% 
% %RMS Error
% for kk = 1:length(u2)
%     temp1 = u2{kk} - u_analyt2{kk};
%     for zz = 1:length(temp1)
%         sums = 0;
%         temp2 = (norm(temp1(zz,:)))^2;
%         sums = sums + temp2;
%     end
%     err2{kk} = sqrt(sums)/length(y0);
% end
% 
% figure; 
% plot(stretches(1:end-1),cell2mat(err2),'o')
% title('RMS Error wrt stretch using kNN Search')
% xlabel('Applied stretch')
% ylabel('RMS')
% 
% %Quiver plots
% figure;
% for jj = 1:length(u2)
%     quiver(stretch_coords(jj).coords(:,1), stretch_coords(jj).coords(:,2),u2{jj}(:,1),u2{jj}(:,2),0)
%     axis image
%     hold on
% end
% xlim([0 1280])
% %%
% 
% % Generate deformed image
% L1 = seedBeadsN(sigma,y1,sizeI);
% 
% %% Save I0 and I1
% 
% save('L0.mat', 'L0');
% save('L1.mat', 'L1');
% 
% %%
% 
% LFM_RayTrace(1280, 800,38,1,.1,15,4,605e-6,500000,y0)