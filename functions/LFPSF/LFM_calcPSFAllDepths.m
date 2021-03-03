% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function psfWaveStack = LFM_calcPSFAllDepths(Camera, Resolution)
%   Computes PSF for all Resolution.depths, exploiting the symetry of Resolution.depths at
%   the same absolute distance to the zero plane

% Offsets the depths in the case of a defocused (2.0) LFM setup; offsetFobj is zero for original (1.0) LFM setup
Resolution.depths = Resolution.depths + Camera.offsetFobj;

psfWaveStack = zeros(length(Resolution.yspace), length(Resolution.xspace), length(Resolution.depths));
disp('Computing PSF for main objective:');
disp('...');
for i = 1:length(Resolution.depths)
    compute_psf = 1;
    idx = 0;
    % Check if the abs(depth) was previoulsy computed, as zero-symetric depths are just conjugates.
    if i > 1
        idx = find(abs(Resolution.depths(1:i-1)) == abs(Resolution.depths(i)));
        if ~isempty(idx)
            compute_psf = 0;
        end
    end
    
    % If depth has not been computed, compute it
    if compute_psf == 1
        tic
        psfWAVE = LFM_calcPSF(0, 0, Resolution.depths(i), Camera, Resolution);
        disp(['PSF: ',num2str(i),'/',num2str(length(Resolution.depths)),' in ',num2str(toc),'s']);
    else
        % if it is exactly the same depth just copy
        if Resolution.depths(i) == Resolution.depths(idx)
            psfWAVE = psfWaveStack(:,:,idx);
        else
            % if it is the negative, conjugate
            psfWAVE = conj(psfWaveStack(:,:,idx));
        end
        disp(['PSF: ',num2str(i),'/',num2str(length(Resolution.depths)),' already computed for depth ',num2str(Resolution.depths(idx))]);
    end
    psfWaveStack(:,:,i)  = psfWAVE;
    
end
disp('...');
end
