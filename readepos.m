function [dx,dy,x,y,z,m,t,vdc,vp,nulls,Nat_pulse]=readepos(file_name)
% [dx,dy,x,y,z,m,t,vdc,vp,nulls,Nat_pulse]=readepos(file_name)
% reads data stored in a epos format called 'file_name.epos' and extracts
% - dx,dy, the detector coordinates
% - x,y,z,m the pos informations
% - t the tof
% - vdc,vp the DC and pulsed voltages
% - nulls the number of pulses between two events
% - Nat_pulse the number of atoms detected on an event
% Author: Baptiste Gault (b.gault@mpie.de)
% tested on MATLAB R2020a

%Authors: Gault B
%(b.gault@mpie.de), 2010
%Australian Centre for Microscopy and Microanalysis (ACMM), The University of Sydney 
%written and tested in MATLAB R2020a 

% BSD 2-Clause License
% 
% Copyright (c) 2010, Gault B
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% opens the file
fid = fopen(file_name, 'r');
disp('File opened...');

%% Reads through the file made of 9 floats separated by 8 (2 imtegers) bytes
lflo=fread(fid, inf, '9*float32',8, 'b');
nb=length(lflo)/9;
%% Makes an array with the list of floats
flo=reshape(lflo, [9 nb]);
%% Creates output
x=flo(1,:);
y=flo(2,:);
z=flo(3,:);
m=flo(4,:);
disp('Good for x,y,z,m...');
t=flo(5,:);
vdc=flo(6,:);
vp=flo(7,:);
disp('Good for tof and voltages...');
dx=flo(8,:);
rdx=rand(1,length(dx));
dx=dx+(rdx-.5)*.05;
dy=flo(9,:);
rdy=rand(1,length(dx));
dy=dy+(rdy-.5)*.05;
disp('Good for detector coordinates..');
clear flo

%% Reinitialise the reading of the file
frewind(fid);
%% Moves after the the first set of floats
fseek(fid,36,-1);
%% Reads through the file made of 2 integers separated by 36 (9 floats) bytes
lin32=fread(fid, inf, '2*int32', 36, 'b');
%% Makes an array with the list of integers
in32=reshape(lin32, [2 nb]);
%% Creates output
nulls=in32(1,:);
Nat_pulse=in32(2,:);
disp('Good for the integers..');
%% Closes the file
fclose(fid);
disp('No worries, file read, and variables created.');

