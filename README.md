# LFM2PT: Light field microscopy multi-particle tracking from high-speed video

LFM2PT is a script bundle that feeds light field microscopy (LFM) data to (2) volumetric particle tracking (PT) to reconstruct 3D volumtric displacements and strains.  

## Purpose
This repository contains the MATLAB m-files for LFM2PT along with synthetic example images. The algorithm relies on two primary precursor code: oLaF (https://gitlab.lrz.de/IP/olaf) and ALT-SCRIPT (early branch from https://github.com/FranckLab/ALT-SCRIPT).

## Running LFM2PT

### Software Required
Development and test occurred in Matlab versions 2019a - 2021a and older versions may work but have not been tested.

The following Toolboxes are required: 
* 'System Identification Toolbox'
* 'Signal Processing Toolbox'
* 'Image Processing Toolbox'
* 'Statistics and Machine Learning Toolbox'
* 'Partial Differential Equation Toolbox'
* 'Wavelet Toolbox'
* 'Curve Fitting Toolbox'
* 'Parallel Computing Toolbox'
* 'MATLAB Parallel Server'
* 'Polyspace Bug Finder'

Workarounds may be possible - attempt at your own risk. 

Other external tools that may be helpful include:
* inpaint_nans and inpaint_nans3d by John D'Errico (2008) (https://www.mathworks.com/matlabcentral/fileexchange?q=inpaint+nans and https://www.mathworks.com/matlabcentral/fileexchange/21214-inpainting-nan-elements-in-3-d)
* shadedErrorBar by Rod Cambell (2009) (https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar and https://github.com/raacampbell/shadedErrorBar)
* sshist by Hideaki Shimazaki (2009) (https://www.mathworks.com/matlabcentral/fileexchange/24913-histogram-binwidth-optimization and http://2000.jukuin.keio.ac.jp/shimazaki)

In addition to several scripts and packages included in the "utils" subfolder.


### Input Images (data)
* 
* 


### Running including example case
1. Make sure that the main files are in the current (working) directory for MATLAB.
2. Copy the desired test images `test_images` directory as needed.
3. Run the `main_LFM2PT.m` script and follow prompts. 


## Files


## FAQ

**

**

## Cite
If used please cite:


## Contact and support
For questions, please first refer to .... Add a new question if similar issue hasn't been reported as a GitHub "Issue". The authors' contact information can be found at [Franck Lab](francklabbackup.me.wisc.edu), via GitHub, or via the associated paper.
