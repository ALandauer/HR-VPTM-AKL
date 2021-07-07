%% Rigid rotation through angle theta about z axis

function [y1] = rotation(y_ref, theta)
%y_ref: original centroid locations in volume [npts x 3] double
%theta: rotation angle in radians


u1(:,1) = y_ref(:,1)*cos(theta) - y_ref(:,2)*sin(theta);
u1(:,2) = y_ref(:,2)*cos(theta) + y_ref(:,1)*sin(theta);
u1(:,3) = y_ref(:,3);

y1 = u1;
end