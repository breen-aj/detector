function [BINS,Edges] = voxelise_2D(x,y, RES)
%this code voxelises up an apt detector hits into 2D bins
%x - x co-ordinates
%y - y-co-ordinates
%RES - the number of bins across the X/Y axis 

%copyright (C) 2022
%author: A.J. Breen, The University of Sydney, 2022
%licence: BSD 2-Clause License - see LICENSE

mi = min(min(x),min(y));
mx = max(max(x),max(y));

Edges{1} = linspace(mi,mx,RES);
Edges{2} = Edges{1};

bin_width = Edges{1}(2) - Edges{1}(1);

BINS = cell(length(Edges{1})-1);

for i = 1: length(x)
    XBin = ceil((x(i) - Edges{1}(1))/bin_width);
    YBin = ceil((y(i) - Edges{1}(1))/bin_width);
    
        %extremes
        if XBin == 0
            XBin = 1;
        end
         if YBin == 0
            YBin = 1;
         end

         if XBin == length(Edges{1})
            XBin = length(Edges{1}) -1;
         end

        if YBin == length(Edges{1})
            YBin = length(Edges{1}) -1;
        end
        
        BINS{XBin, YBin} = [BINS{XBin, YBin} i];

end

end
    

