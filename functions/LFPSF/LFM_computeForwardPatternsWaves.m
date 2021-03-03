% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function H = LFM_computeForwardPatternsWaves(psfWaveStack, MLARRAY, Camera, Resolution)

% ComputeForwardPatternsWaves: Compute the forward projection for every source point (aa,bb,c) of our square shaped patch around the central microlens
% Take the PSF incident on the MLA (psfWAVE_STACK), pass it through the
% microlens array (MLARRAY), and finally propagate it to the sensor

% for regular grids compute the psf for only one quarter of coordinates (due to symmetry)
if strcmp(Camera.range, 'quarter')
    coordsRange  = Resolution.TexNnum_half;
else
    coordsRange  = Resolution.TexNnum;
end

Nnum_half_coord = Resolution.TexNnum_half ./ Resolution.texScaleFactor;
sensorRes = Resolution.sensorRes;

% Resolution.texScaleFactor(1/2) is actually (Resolution.texRes(1/2) * M / Resolution.sensorRes(1/2))^-1  
H = cell( coordsRange(1), coordsRange(2), length(Resolution.depths) ); 
for c = 1:length(Resolution.depths)
    psfREF = psfWaveStack(:,:,c);
    for aa_tex = 1:coordsRange(1)
            aa_sensor = aa_tex / Resolution.texScaleFactor(1);
        parfor bb_tex = 1:coordsRange(2)
            bb_sensor = bb_tex / Resolution.texScaleFactor(2);
            
            % shift the native plane PSF at every (aa, bb) position (native plane PSF is shift invariant)
            psfSHIFT = imShift2(psfREF, round((aa_sensor-Nnum_half_coord(1))), round((bb_sensor-Nnum_half_coord(2))) );
            
            % MLA transmittance
            psfMLA = psfSHIFT.*MLARRAY;          
            
            % propagate the response to the sensor via Rayleigh-Sommerfeld diffraction
            LFpsfAtSensor = prop2Sensor(psfMLA, sensorRes, Camera.mla2sensor, Camera.WaveLength, 0);
            
            % shift the response back to center (we need the patterns centered for convolution)
            LFpsf = imShift2(LFpsfAtSensor, round(-(aa_sensor-Nnum_half_coord(1))), round(-(bb_sensor-Nnum_half_coord(2))) );
            
            % store the response pattern 
            H{aa_tex,bb_tex,c} = sparse(abs(double(LFpsf).^2));
        end
    end
    disp(['Forward Patterns, depth: ', num2str(c), '/', num2str(length(Resolution.depths))]);
end