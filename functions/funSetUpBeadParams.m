function [beadParameter] = funSetUpBeadParams(beadParameter)

% Define default values
thres = 0.5;
minSize = 5;
maxSize = Inf;
winSize = [15,15,15];

beadSize = 3;          % Estimated radius of a single particle
numBeadsPSF = 1;       % number of bead to use for mean psf (method 3)
PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
distMissing = 5;       % Distance threshold to check whether particle has a match or not 
color = 'white';       % By default
detectionMethod = 1;    % Default
saveIntermediates = 1;

dccd = [1,1,1];
abc = [1,1,1];
forloop = 1;
randNoise = 1/10^7; % Something small
xy = 1; % Assume symmetrical if not given
z = 1;  % Assume symmetrical if not given

diameter = 5;   % Assume if not given
ratThresh = 0.2;
circThresh = 1.0;
smoothFac = 0.15;

%deconv params
fileFolder = './';
deconvThresh = 0.05;
deconvPrefilter = true; %true/false gaussian prefilter option
deconvIter = 5;
psfSize = [40,40]; %[x,y] size, z is picked by user

p = inputParser;
addParameter(p,'thres',thres);
addParameter(p,'minSize',minSize);
addParameter(p,'maxSize',maxSize);
addParameter(p,'winSize',winSize);
addParameter(p,'beadSize',beadSize);
addParameter(p,'numBeadsPSF',numBeadsPSF);
addParameter(p,'PSF',PSF);
addParameter(p,'distMissing',distMissing);
addParameter(p,'color',color);
addParameter(p,'detectionMethod',detectionMethod);
addParameter(p,'dccd',dccd);
addParameter(p,'abc',abc);
addParameter(p,'forloop',forloop);
addParameter(p,'randNoise',randNoise);
addParameter(p,'xy',xy);
addParameter(p,'z',z);
addParameter(p,'diameter',diameter);
addParameter(p,'ratThresh',ratThresh);
addParameter(p,'circThresh',circThresh);
addParameter(p,'smoothFac',smoothFac);
addParameter(p,'deconvThresh',deconvThresh);
addParameter(p,'deconvPrefilter',deconvPrefilter);
addParameter(p,'deconvIter',deconvIter);
addParameter(p,'psfSize',psfSize);
addParameter(p,'fileFolder',fileFolder);
addParameter(p,'saveIntermediates',saveIntermediates);

parse(p,beadParameter)

beadParameter = p.Results;

