This repository contains the custom MATLAB code used in the study: Extended wakefulness fragments memory-consolidating sleep in aging.

Software: MATLAB R2024b (MathWorks) following the official installation instructions provided by MathWorks.

Operating Systems: The code was developed and tested using MATLAB R2024b on a standard desktop computer.

Additional Dependencies: The code requires only standard MATLAB R2024b functions and toolboxes. All custom functions developed for this study are included in this repository. No additional third-party software, custom packages, or external dependencies are required.
To load fiber photometry data we use a custom function provided by Tucker Davis Technologies called "TDTbin2mat" that you can find in their repository (https://github.com/tdtneuro/TDTMatlabSDK). To load EEG and EMG data obtained using Sleepscore (ViewPoint, France) we use a custom function provided by ViewPoint Behavior Technology called "loadEXP" - this function can be obtained upon request (www.viewpoint.fr). The infraslow power analysis is adpated from ISr morlet transform function (https://github.com/arnaudboutin/Spindle-SO-package).

The code is organized into folders with "scripts" containing scripts used for analysis and "functions" containing custom functions that are used in the scripts. 

Demo Dataset: The dataset contains representative EEG, NE, and full hypnogram data that can be used to test the analysis pipeline and verify successful installation. For the testing analysis, you need to load mat file on matlab, add the repository folder and all subfolders to the MATLAB path, run the analysis script according to the instructions provided in the repository.
