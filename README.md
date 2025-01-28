# detector
This repository contains matlab code for sophisticated multimetric APT detector mapping. 

Please refer to the following journal article for details:

Breen AJ, Day AC, Lim B, Davids WJ, Ringer SP. Revealing latent pole and zone line information in atom probe detector maps using crystallographically correlated metrics. Ultramicroscopy. 2023 Jan 1;243:113640.

[https://doi.org/10.1016/j.ultramic.2022.113640](https://doi.org/10.1016/j.ultramic.2022.113640)

Please cite this article if you use the code for any presented/published works.

# How to use the code

Download the repository and add the code to your read path in MATLAB.

## Mapping density, electric field distribution, distance between successive evaporation events and multiple hits 
Open the script: 'APT_DETECT_GITHUB.m' in the editor. Have the .epos file of interest in the current folder/working directory. Change the file name in line 10 of APT_DETECT_GITHUB.m to your .epos file of interest. Change the detection event sequence start and end points in line 35 and 36. Usually a sequence of 1-2 million ions from a single grain works best. However, for the electric field distribution map, a larger number of sequential detection events for each grain of interest (~ > 10 million) performs better. For nanocrystalline grains you may want to use less than 1 million sequential ions to ensure all ions are coming from the same grain - otherwise there will be overlap in the crystallographic patterns. There are also numerous points throughout this main script (and the called functions) where you may want to modify the code. Please follow the in-script and in-function comments for additional information.

The code is currently configured to perform detector mapping on the listed metrics above for 2 million sequential ion hits. Test data used in the journal article listed above can be found at:

[test data](https://cloudstor.aarnet.edu.au/plus/s/uQ59EGNXxxv2g2C)

## Spatial signal mapping
Open the main script "SPATIAL_SIGNAL_DF_MAP.m". Have the .epos file of interest in the current folder. Change the file name in line 29 to your .epos file of interest. Change lines 40-41 to set the sequence range of interest (typically 1-2 million ions). Follow the in-script and in-function comments for additional modifications and changes.

The script is configured to do a spatial signal map on the Al atoms contained within 1-2 million sequential detection events from the gamma prime phase of the example dataset shared above. 

# Other information
The code has been built and tested on MATLAB R2020a.
The code can serve as a template for developing crystallographic spatial signal maps across the detector space for APT data. 
For more information/help please email andrew.breen@sydney.edu.au
