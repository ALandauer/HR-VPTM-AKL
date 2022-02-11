function [vol,BeadPara] = funPreprocLocalizeAC(vol_in,BeadPara,file_names,ImgSeqNum)
%Function to gather a PSF, run deconv, and collect threshold info
%
%--- INPUTS ---
%  vol_in      : input volume to deconv
%  BeadPara    : bead localization parameters struct
%  ImgSeqNum   : volume number in the sequence
%
%--- OUTPUTS ---
% vol: final deconvolved reconstruction
% thresh_est: threshold estimate from the deconv'd image
%
%Alex Landauer, 2020-10-14
%Franck Lab, Brown Univerisity, University of Wisc - Madison; NIST MML

deconv_thresh = BeadPara.deconvThresh;
deconv_prefilter = BeadPara.deconvPrefilter;
deconv_iter = BeadPara.deconvIter;
fileFolder = BeadPara.fileFolder;
cur_deconv_img = [fileFolder,filesep,'deconv_imgs',filesep,'deconv_',file_names{ImgSeqNum}];

if ~exist([fileFolder,filesep,'deconv_imgs'],'dir')
    mkdir([fileFolder,filesep,'deconv_imgs'])
end

if BeadPara.saveIntermediates == 1
    if exist(cur_deconv_img,'file') == 2
        vol_ = load(cur_deconv_img);
        vol=vol_.vol;
    else
        if exist([fileFolder,filesep,'PSF.mat'],'file') == 2
            load([fileFolder,filesep,'PSF.mat'])
            BeadPara.PSF = PSF;
        else
            PSF = funGetPsfMultibead(vol_in,BeadPara.psfSize,BeadPara.numBeadsPSF);
            save([fileFolder,filesep,'PSF.mat'],'PSF')
            BeadPara.PSF = PSF;
        end
        
        %prefilter with a small Gaussian
        if deconv_prefilter
            vol_in = imgaussfilt3(double(vol_in));
            psf_in = imgaussfilt3(double(PSF));
        else
            psf_in = PSF;
        end
        
        %normalize input volume and psf to range [0,1]
        vol_deconv = vol_in/max(vol_in(:));
        vol_deconv(vol_deconv < deconv_thresh) = 0;
        
        psf = psf_in/max(psf_in(:));
        psf(psf<deconv_thresh) = 0;
        
        %%%%%%%%%%%%%for testing
        % psf = imnoise(psf);
        
        %Do the deconvolution
        disp('%%%%%% Starting deconvolution %%%%%%')
        
        %L-R type deconv
        vol = deconvlucy(vol_deconv, psf, deconv_iter);
        
        %renormalize output
        vol = vol/max(vol(:));
        vol(vol < deconv_thresh) = 0;
        
        save(cur_deconv_img,'vol');
        
    end
    
    
else
    if exist([fileFolder,filesep,'PSF.mat'],'file') == 2
        load([fileFolder,filesep,'PSF.mat'])
        BeadPara.PSF = PSF;
    else
        PSF = funGetPsfMultibead(vol_in,BeadPara.psfSize,BeadPara.numBeadsPSF);
        save([fileFolder,filesep,'PSF.mat'],'PSF')
        BeadPara.PSF = PSF;
    end
    
    %prefilter with a small Gaussian
    if deconv_prefilter
        vol_in = imgaussfilt3(double(vol_in));
        psf_in = imgaussfilt3(double(PSF));
    else
        psf_in = PSF;
    end
    
    %normalize input volume and psf to range [0,1]
    vol_deconv = vol_in/max(vol_in(:));
    vol_deconv(vol_deconv < deconv_thresh) = 0;
    
    psf = psf_in/max(psf_in(:));
    psf(psf<deconv_thresh) = 0;
    
    %%%%%%%%%%%%%for testing
    % psf = imnoise(psf);
    
    %Do the deconvolution
    disp('%%%%%% Starting deconvolution %%%%%%')
    
    %L-R type deconv
    vol = deconvlucy(vol_deconv, psf, deconv_iter);
    
    %renormalize output
    vol = vol/max(vol(:));
    vol(vol < deconv_thresh) = 0;
    
    
end

%Matlab 'blind' type deconv
% V = 0.001; %amount of noise supression
% wt = zeros(size(vol_deconv));
% wt(10:end-9,10:end-9) = 1; %weighting function for px used

% [vol, psfr] = deconvblind(vol_deconv,psfi,iter,V,wt);

% vol = deconvreg(I, PSF); %a third possible deconv

% get the threshold estimate if this is the first image in the sequence
if ImgSeqNum == 1
    f1 = figure;
    histogram(vol(vol~=0))
    drawnow
    f2 = figure;
    imagesc3D(vol),colorbar
    drawnow
    
    BeadPara.thres = input('Enter binarization threshold estimate: ');
    try
        close(f1)
        close(f2)
    catch
    end
    YN = 0;
    while YN == 0
        f3 = figure;
        imagesc3D(vol>BeadPara.thres),colorbar
        drawnow
        
        YN_ = input('Binarization okay? (Y/N) [N]: ','s');
        
        if strcmpi(YN_,'Y')
            YN = 1;
        else
            f1 = figure;
            histogram(vol(vol~=0))
            drawnow
            f2 = figure;
            imagesc3D(vol),colorbar
            drawnow
            BeadPara.thres = input('Enter binarization threshold estimate from histogram: ');
        end
        
    end
    try
        close(f1)
    catch
    end
    try
        close(f2)
    catch
    end
    try
        close(f3)
    catch
    end
    %         thres_est = 0.075;
end

disp('%%%%%% Deconvolution complete! %%%%%%')


