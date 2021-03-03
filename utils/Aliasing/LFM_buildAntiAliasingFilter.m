% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function lanczos2FFT = LFM_buildAntiAliasingFilter(filterSize, widths, n)

[x, y] = meshgrid(-floor(filterSize(2)/2):floor(filterSize(2)/2),...
        -floor(filterSize(1)/2):floor(filterSize(1)/2)); %% note the [x,y] convention

% compute the Lanczos windowed sinc kernel 
kernelSinc = zeros(filterSize);
lanczos2 = zeros(filterSize);
lanczos2FFT = zeros(filterSize);
for i = 1:size(kernelSinc,3)
    
    % ideal filter
    fy = 1/widths(i,1);
    fx = 1/widths(i,2);
    x_f = x*fx;
    y_f = y*fy;
    kernelSinc(:,:,i) = (sinc(x_f).*sinc(y_f));
    
    % nX Lanczos window -- depending on the texture features
    current_window = sinc(x_f/n).*sinc(y_f/n); % ones(size(current_kernel)); 
    current_window(sqrt(y.^2 + x.^2) > 3*widths(i,1)) = 0;
    
    % final windowed ideal filter
    lanczos2(:,:,i) = kernelSinc(:,:,i).*current_window;
    current_kernel = lanczos2(:,:,i);
    lanczos2(:,:,i) = (lanczos2(:,:,i)./sum(current_kernel(:)));
    lanczos2FFT(:,:,i) = fft2(lanczos2(:,:,i));
end