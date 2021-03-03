%% Compute error

load(exactFileName)

parCoordA_exact = x0{1};
parCoordB_exact = x1{1};
u_exact = u{1};

f_o_s = 3; 
[matchesA] = f_track_nearest_neighbour3( parCoordA, parCoordA_exact, f_o_s );
[matchesB] = f_track_nearest_neighbour3( parCoordB, parCoordB_exact, f_o_s );

% Assemble variable "track"
trackA = zeros(size(parCoordA,1),1);
for tempi = 1:size(matchesA)
   trackA(matchesA(tempi,1)) = matchesA(tempi,2); 
end

trackB = zeros(size(parCoordB,1),1);
for tempi = 1:size(matchesB)
   trackB(matchesB(tempi,1)) = matchesB(tempi,2); 
end


%% trackA2B_exact
track_A2B_exact = zeros(size(parCoordA,1),1);
for tempi = 1:size(parCoordA,1)
   [row,~] = find ( trackB == trackA(tempi) );
    try 
        track_A2B_exact(tempi) = row; 
    catch
    end
    
end

[row1,~] = find(track_A2B>0);
[row,~] = find(abs(track_A2B(row1)-track_A2B_exact(row1))>0);
% length(row1)
% length(row)
MismatchRatio = length(row)/length(row1);
disp( ['Mismatch ratio: ',  num2str( length(row) ),' / ', num2str(length(row1)), ' = ', num2str(MismatchRatio) ]);
    

%% u_A2B_exact
u_A2B_exact = u_exact(  trackA(track_A2B>0), :);

parCoordA_ = parCoordA( (track_A2B>0),:);  

% %%%%% Exact disp plot in reference configuration %%%%%
% figure, quiver3(parCoordA_(:,1),parCoordA_(:,2),parCoordA_(:,3),u_A2B_exact(:,1),u_A2B_exact(:,2),u_A2B_exact(:,3));
% set(gca,'fontsize',30); view(2); box on; axis equal; axis tight; % view(3);
% title('Exact disp field','fontweight','normal');

% %%%%% Exact disp plot in deformed configuration %%%%%
figure, quiver3(parCoordA_(:,1)+u_A2B_exact(:,1), parCoordA_(:,2)+u_A2B_exact(:,2), ...
    parCoordA_(:,3)+u_A2B_exact(:,3), u_A2B_exact(:,1),u_A2B_exact(:,2),u_A2B_exact(:,3));
set(gca,'fontsize',30); view(2); box on; axis equal; axis tight; % view(3);
title('Exact disp field','fontweight','normal');

clear u_A2B_err
if iterNum == 1
    u_A2B_err = u_A2B - u_A2B_exact;
else
    u_A2B_err = u_A2B - u_B2A_curr(track_A2B(track_A2B>0),:) - u_A2B_exact;
end
 
% %%%%% Error disp plot in reference configuration %%%%%
% figure, quiver3(parCoordA_(:,1),parCoordA_(:,2),parCoordA_(:,3),u_A2B_err(:,1),u_A2B_err(:,2),u_A2B_err(:,3));
% set(gca,'fontsize',30); view(2); box on; axis equal; axis tight; % view(3);

% ----- histogram -----
u_A2B_err = sqrt( sum( u_A2B_err.^2, 2) );

figure, h = histogram( u_A2B_err,100 ); set(gca,'fontsize',30);

