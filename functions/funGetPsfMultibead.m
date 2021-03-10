function [psf_mean] = funGetPsfMultibead(img,psf_size,num_beads)
% Function to gather a PSF of a bead from a post-reconstruction image
%
% --- INPUTS ---
%  vol_in    : input volume image
%  psf_size  : approximate (due to rounding) size of the output
%                     point spread function (PSF)
%  num_beads : number of beads to segment to get mean psf
%
% --- OUTPUTS ---
% psf: single bead image from the reconstruction
%
% Alex Landauer, 2020-10-14
% Franck Lab, Brown Univerisity, University of Wisc - Madison; NIST MML

%normalize the image to max = 1
img = double(img)/max(double(img(:)));

disp('%%%%%% Collect bead PSF from volume %%%%%%')

%display for the user and prompt for input
f1 = figure;
try
    imagesc3D(img)
catch
    imshow3D(img)
end
drawnow

%get z-plane for each
for ii = 1:num_beads
    z_bottom_plane(ii) = input('Enter bead bottom plane: ');
    z_top_plane(ii) = input('Enter bead top plane: ');
    z_plane(ii) = z_bottom_plane(ii) + round((z_top_plane(ii) - z_bottom_plane(ii))/2);
end
try
    close(f1)
catch
end

disp('Click on bead centers in order')
for ii = 1:num_beads
    
    %select x,y locations
    h = figure;
    z_slice = img(:,:,z_plane(ii))/max(max(img(:,:,z_plane(ii))));
    imagesc(z_slice),axis image
    drawnow
    
    [y_center(ii),x_center(ii)] = ginput(1);
    close(h)
    
    %get the range(s) within which the bead(s) resides
    range_x(:,ii) = (round(x_center(ii)-psf_size(1)/2)):(round(x_center(ii)+psf_size(1)/2));
    range_y(:,ii) = (round(y_center(ii)-psf_size(2)/2)):(round(y_center(ii)+psf_size(2)/2));
    range_z{ii} = z_bottom_plane(ii):z_top_plane(ii);
    
    %collect the psf
    psfs_{ii} = double(img(range_x(:,ii),range_y(:,ii),range_z{ii}));
    
end

for ii = 1:num_beads
    z_b_dif = z_bottom_plane(ii) - min(z_bottom_plane);
    z_t_dif = max(z_top_plane) - z_top_plane(ii);
    psfs_pad{ii} = cat(3,nan*ones(length(range_x),length(range_y),z_b_dif),...
        psfs_{ii}(:,:,:),nan*ones(length(range_x),length(range_y),z_t_dif));
    psfs(:,:,:,ii) = psfs_pad{ii}(:,:,:);
end

%find mean
psf_mean = nanmean(psfs,4);

%display the result
h = figure;
imagesc3d(psf_mean)
drawnow

%verify with user if it is okay
prompt = 'PSF ok? Y/N [Y]: ';
yn = input(prompt,'s');
if isempty(yn)
    yn = 'Y';
end

if yn == 'y' || yn == 'Y'
    try
    close(h)
    catch
    end
else
    error('psf not accepted')
end

try
    close(f1)
catch
end

















