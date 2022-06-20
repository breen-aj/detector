function [x, y, z, R_start, R_end] = atomProbeRecon05(detx, dety,h,k, V, kf, ICF, avgDens, Fevap, flightLength, detEff)
%atom probe reconstruction after: Gault et al., Ultramicroscopy 111 (2011) 448 - 457
% detx, dety are the detector hit coordinates in mm
% h and k are the x and y centres of the detector ROI
% V is the volage (standing voltage + voltage pulse) that the ions were
% evaporated at
% kf is the field factor and ICF is the image compression factor
% ICF is the image compression factor
% avgDens is the theoretical average density of the specimen being analysed
% in V/nm^3
% Fevap is the average evaporation field of the ions in V/nm
% flightLength is the flight path length that the APT data was collected at
% detEff is the detector efficiency of the APT instrument used to collect
% the data

% author: Peter Felfer, Australian Centre for Microscopy and Microanalysis,
% The University of Sydney, ~ 2015

% with minor modifications by Andrew Breen, Australian Centre for Microscopy and Microanalysis, The University of Sydney, 2022

% licence: BSD 2-Clause License - see LICENCE

%% constants and variable setup

% detector coordinates in polar form
[ang, rad] = cart2pol(detx-h, dety-k);

%max(rad)

% calcualting effective detector area:
Adet = (max(rad))^2 * pi();

% radius evolution from voltage curve (in nm)
Rspec = V/(kf * Fevap);

%output start and end radius
R_start = Rspec(1);
R_end = Rspec(end);

%% calcualte x and y coordinates

% compressed angle
thetaP = atan(rad / flightLength); % mm/mm

% launch angle
theta = thetaP + asin((ICF - 1) * sin(thetaP));

% distance from axis and z shift of each hit
[zP, d] = pol2cart(theta, Rspec); % nm

% x and y coordinates from the angle on the detector and the distance to
% the specimen axis.
[x, y] = pol2cart(ang, d); % nm



%% calculate z coordinate
% the z shift with respect to the top of the cap is Rspec - zP
zP = Rspec - zP;

% accumulative part of z
omega = 1 ./ avgDens; % average atomic volume of each atom in nm^3

dz = omega * flightLength^2 * kf^2 * Fevap^2 / (detEff * Adet * ICF^2) * V.^-2; % nm^3 * mm^2 * V^2/nm^2 / (mm^2 * V^2) 

% wide angle correction
z = cumsum(dz) + zP;

%% end of reconstruction
end


