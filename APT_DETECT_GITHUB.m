%this script combines correlative crystallography metrics of the detector hit maps

%copyright (C) 2022
%author: Andrew Breen, The University of Sydney, 2022
%licence: BSD 2-Clause License - see LICENSE

close all
clear all
clc



%read in the APT dataset
[dx,dy,x,y,z,m,t,vdc,vp,nulls,Nat_pulse]=readepos('R18_59377-v01.epos'); %AM In738 test file. Please replace with your own .epos file 
 
%% clip file

rad = 40; %crop the detector map to a desired radius

keep = dx.^2 + dy.^2 < rad.^2;

dx2 = dx(keep);
dy2 = dy(keep);
x2 = x(keep);
y2 = y(keep);
z2 = z(keep);
m2 = m(keep);
t2 = t(keep);
vdc2 = vdc(keep);
vp2 = vp(keep);
nulls2 = nulls(keep);
Nat_pulse2 = Nat_pulse(keep);

% Change the ion sequence range. Usually 1-2 million sequential atoms is ideal.
% Exceptions are nanocrystalline datasets (use a small range to capture a
% single grain
%electric field map: this works better when a large number of ions are
%considered i.e. >> 10 million
start_num = 4e6; 
end_num = 6e6;


dx3 = dx2(start_num:end_num);
dy3 = dy2(start_num:end_num);
x3 = x2(start_num:end_num);
y3 = y2(start_num:end_num);
z3 = z2(start_num:end_num);
m3 = m2(start_num:end_num);
t3 = t2(start_num:end_num);
vdc3 = vdc2(start_num:end_num);
vp3 = vp2(start_num:end_num);
nulls3 = nulls2(start_num:end_num);
Nat_pulse3 = Nat_pulse2(start_num:end_num); 

% %% run maps

sm = 0; %smooth
%% this calculates density (1)
[N1,a1,b1] = DENSITY_MAP(dx3,dy3,m3, sm, 0);

%% this calculates Nulls (2)
[N2, a2, b2] = NULL_MAP(dx3, dy3, nulls3, sm);

%% this calculates multiple ion hits (3)
[N3, a3, b3] = MULTI_MAP(dx3, dy3, Nat_pulse3, sm);

%% this calculates density of a particular species
%ranges Al species R18_59377. Update for your ion charge states
RANGES = [8.847, 9.146; 13.372, 13.898]; %Al 3+, 2+
%RANGES = [13.347, 13.898; 26.855, 27.357]; %Al 2+, 1+


%% this calculated electric field map based on RANGES
[N5,a5,b5] = EFIELD_MAP(dx2, dy2, m2, sm, RANGES);

%% this calculates distance between evaporation events
% steps is number of atoms to measure distance to, either side of evaporated atom 
% more steps - good for ring counting
% less stesp - good for pole clarity
%typically choose between 2 - 20
steps = 2;

[N7,a7,b7] = EVAPSEQ_MAP(dx3, dy3, m3, sm, 0, steps);