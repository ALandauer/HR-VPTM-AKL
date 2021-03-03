% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function newMask = LFM_fixMask(mask, NewLensletSpacing, gridType)
% This function corrects a mask for holes and overlaps
% Uncomment plot related lines for a visual hint

trial_space = zeros(3*size(mask));
r_center = ceil(size(trial_space,1)/2);
c_center = ceil(size(trial_space,2)/2);
rs = NewLensletSpacing(1);
cs = NewLensletSpacing(2);

if (strcmp(gridType, 'hex'))
    trial_space([r_center - rs, r_center + rs],[c_center - round(cs/2), c_center + round(cs/2)]) = 1;
    trial_space([r_center, r_center, r_center],[c_center - cs, c_center, c_center + cs]) = 1;
else
    if (strcmp(gridType, 'reg'))
        trial_space([r_center - rs, r_center + rs, r_center],[c_center - cs, c_center + cs, c_center]) = 1;
    end
end

space = conv2(trial_space, mask, 'same');
[r, c] = size(mask);
space_center = space(r+1:r*2, c+1:c*2);

% figure;
newMask = mask; 
for i = 1:r
    for j = 1:c
        % fix holes
        if(space_center(i,j) == 0)
            newMask(i,j) = 1;
            space = conv2(trial_space, newMask, 'same');
            space_center = space(r+1:r*2, c+1:c*2);
           
%             subplot(2,2,3); imagesc(space_center); title('central part of the space above')
%             subplot(2,2,1:2);imagesc(space); title('3x3 test space')
%             subplot(2,2,4);imagesc(new_mask); title('new mask')
%             pause
        end
        
        % fix overlap
        if(space_center(i,j) == 2)
            newMask(i,j) = 0;
            space = conv2(trial_space, newMask, 'same');
            space_center = space(r+1:r*2, c+1:c*2);
%              
%             subplot(2,2,3); imagesc(space_center); title('central part of the space above')
%             subplot(2,2,1:2);imagesc(space); title('3x3 test space')
%             subplot(2,2,4);imagesc(new_mask); title('new mask')
%             pause
        end
        
    end
end