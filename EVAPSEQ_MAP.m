function [N_new,a,b] = EVAPSEQ_MAP(dx, dy, m, sm, RANGES, steps)
%This function plots a field desorption map of inputted detector
%co-ordinates. Note slices ~ 1-2 million atoms usually show clearest FDM
%parent function: correlative detector

%dx - inputted X co-ordinates
%dy - inputted Y co-ordinates
%sm - smooth SDM? (1 = yes, 0 = no)
%RANGES -filter dx, dy on inputted - if 0 skip
%Steps - the numer of evaporated ion steps to measure distance to. 

%copyright (C) 2022
%author: A.J. Breen, The University of Sydney, 2022
%licence: BSD 2-Clause License - see LICENCE


fsize = 9;

%% filter dx and dy based on chemistry
 dx02 = [];
 dy02 = [];

if RANGES ~= 0
    range_num = size(RANGES, 1);
    for i = 1:range_num
        dx02 = [dx02 dx(m > RANGES(i, 1) & m < RANGES(i, 2))];
        dy02 = [dy02 dy(m > RANGES(i, 1) & m < RANGES(i, 2))];
    end
else
    dx02 =dx;
    dy02 =dy;
    
end
        
%% calculate the distance between evaporation steps

dist_sq = zeros(length(dx02)-steps*2, steps*2+1);

for i = 1: length(dx02)-2*(1+steps)
    for j = 1:steps*2+1
        dist_sq(i,j) = (dx(i+steps) - dx(i+j-1)).^2 + (dy(i+steps) - dy(i+j-1)).^2; %calcuating dist squared is faster than dist - don't use sqrt
    end
end

%remove atoms counting themselves
dist_sq(:,steps+1) = [];
        
dist_ave = sqrt(sum(dist_sq, 2)./size(dist_sq,2));%dist average

%voxelise and calculate the score of each bin based on the ave dist to next
%evaporated ions in that bin. 

%remove atoms on edges of array

dx = dx(steps+1:length(dx)-steps);
dy = dy(steps+1:length(dy)-steps);

mi = min(min(dx),min(dy));
mx = max(max(dx),max(dy));

%using rices_rule
%RES = floor(2*length(dx)^(1/3));

%override
RES = 100;

[BINS,Edges] = voxelise_2D(dx,dy, RES);

Ctr{1} = Edges{1} + (Edges{1}(2) - Edges{1}(1))/2;
Ctr{1}(end) = [];
Ctr{2} = Ctr{1};

[b,a] = meshgrid(Ctr{1}, Ctr{2});

%find average score for each voxel

N = zeros(size(b,1), size(b,2));

for i = 1:size(b, 1)
    for j = 1:size(b,2)
        if isempty(BINS{i,j}) == 0 
            N(i,j) = sum(dist_ave(BINS{i,j}))/length(BINS{i,j});
        else
            N(i,j) = 0;
        end
    end
end


if sm == 1
    N = smoothn(N, 0.1);
end 


%% remove background hump

% normalise the plots and combine

Z = reshape(N,[],1); %z co-ord
X = reshape(a,[],1); %x co-ord
Y = reshape(b,[],1); %y-co-ord

X = X(Z~=0);
Y = Y(Z~=0);
Z = Z(Z~=0);

%note: only fit to nice atoms in the centre e.g. 10-20 mm out to avoid edge
%effects
rad = 25;

X2 = X(X.^2 +Y.^2 < rad.^2);
Y2 = Y(X.^2 +Y.^2 < rad.^2);
Z2 = Z(X.^2 +Y.^2 < rad.^2);


sf = fit([X2, Y2], Z2, 'poly22');
N_hump = sf.p00 + sf.p10*a + sf.p01*b + sf.p20*a.^2 + sf.p11*a.*b + sf.p02*b.^2; 

N_new = N-N_hump;
N_new(N == 0) = 0;
N_new(N == 0) = min(min(N_new));
N_new = N_new-min(min(N_new));
%N_new = N_new.^2;

%manual removal of outliers
%histogram the distribution,spread as little as possible from 0
% figure
% histogram(N_new, -1000:10:1000);
% %disp('break here');
% N_new(N_new > 220) = 0;
% N_new = N_new+850;
% N_new(N_new<0) = 0;
% N_new(N==0) = 0;


%remove outliers
%hist(N_new)

%% plot figure

fig_dim_x = 6;%change the dimensions of the figure here
fig_dim_y = 5;

fig = figure('Units', 'centimeters');
fig.Position = [5,5,fig_dim_x, fig_dim_y];

h = surf(a,b, N_new);

%colormap gray
colormap inferno
shading interp
view(2)
axis square
grid off
cb = colorbar;
set(get(cb,'title'), 'string', 'S (mm)', 'fontsize', 9);
cb.FontSize = 9;
ax1 = gca;
set(ax1,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on')
axis off
xlim([min(dx) max(dx)])
ylim([min(dy) max(dy)])

set(gcf, 'PaperPositionMode','auto');
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'color', 'white');
%axis off

%change caxis

mid_c = mean(N_new(a.^2 +b.^2 < rad.^2));

std_c = std(N_new(a.^2 +b.^2 < rad.^2));

caxis([mid_c - 2*std_c, mid_c + 2*std_c]); 

end

