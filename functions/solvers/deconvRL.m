% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu

function recon = deconvRL(forwardFUN, backwardFUN, img, iter, init)
% Richardson-Lucy deconvolution algorithm
fprintf('\nDeconvolution:')

% Initialize volume
recon = init; 

for i = 1:iter
    tic
    % simulate forward projection of the current reconstructed volume 
    fpj = forwardFUN(recon);
    
    % compute error towards the real image
    errorBack = img./fpj;
    
    % make sure the computations are safe
    errorBack(isnan(errorBack)) = 0;
    errorBack(isinf(errorBack)) = 0;
    errorBack(errorBack < 0) = 0;
    
    % backproject the error
    bpjError = backwardFUN(errorBack);

    % update the result
    recon = recon.*bpjError;
    ttime = toc;
    fprintf(['\niter ' num2str(i) ' | ' num2str(iter) ', took ' num2str(ttime) ' secs']);
end