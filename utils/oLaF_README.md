# oLaF - A flexible 3D reconstruction framework for light field microscopy

Copyright (c)2017-2020 Anca Stefanoiu, Josue Page, Tobias Lasser

## Overview
The light field microscope (LFM) enables scan-less 3D imaging of fluorescent specimens by incorporating an array of micro-lenses (MLA) into the optical path of a conventional wide-field microscope. Thus, both spatial and directional light field information is captured in a single shot, allowing for subsequent volumetric reconstruction of the imaged sample.

oLaF is a flexible and efficient Matlab framework for 3D reconstruction of light field microscopy data. It is designed to cope with various LFM configurations in terms of MLA type (regular vs. hexagonal grid, single-focus vs. mixed multi-focus lenslets) and placement in the optical path (original 1.0 vs. defocused 2.0 LFM vs. Fourier LFM designs).

## Acknowledgements
oLaF evolved gradually improving and expanding on the functionality in [1]. Where possible, the naming conventions were kept for the sake of relatability and convenience of the shared users.      
The pre-processing of the raw light field images relies on the image rectification functionality [2] in the Light Field Toolbox for Matlab by Donald G. Dansereau.

*[1] Supplementary software: R. Prevedel et. al., "Simultaneous whole-animal 3D-imaging of neuronal activity using light field microscopy", Nat. Methods 11, 727–730 (2014).*

*[2] D. G. Dansereau, O. Pizarro, and S. B. Williams, "Decoding, calibration and rectification for lenselet-based plenoptic cameras", Computer Vision and Pattern Recognition (2013).*

## Getting started
Code/mainLFM.m serves as a demo script of the conventional LFM functionality in oLaF v3.0.

Code/mainFLFM.m serves as a demo script of the Fourier LFM functionality in oLaF v3.0. 

SampleData/ folder contains various sample datasets. 

## Documentation
https://mediatum.ub.tum.de/doc/1522002/file.pdf

## Citing
When using oLaF, please reference the following citations:

*[1] A. Stefanoiu et. al., "Artifact-free deconvolution in light field microscopy", Opt. Express, 27(22):31644, (2019).*

*[2] A. Stefanoiu et. al., "What about computational super-resolution in fluorescence Fourier light field microscopy?", Opt. Express, 28, 16554 (2020).*

*[3] A. Stefanoiu et. al., "Deconvolution in Fourier integral microscopy", Proc. SPIE 11396, Computational Imaging V (2020).*

## Feedback

Feedback and suggestions are always welcome. 

anca.stefanoiu *(at)* tum *(dot)* de / josue.page *(at)* tum *(dot)* de.
