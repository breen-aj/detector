function [N,a,b] = DENSITY_MAP(dx, dy, m, sm, RANGES)
%This function plots a field desorption map of inputted detector
%co-ordinates. Note slices ~ 1-2 million atoms usually show clearest FDM
%parent fucntion: correlative detector

%dx - inputted X co-ordinates
%dy - inputted Y co-ordinates
%m - inputted mass/charge of ions
%sm - smooth SDM? (1 = yes, 0 = no)
%RANGES -filter dx, dy on inputted - if 0 skip

%copyright (C) 2022
%author: A.J. Breen, The University of Sydney, 2022
%licence: BSD 2-Clause License - see LICENSE


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
        
    

%% calculating bin centers for FDM


mi = min(min(dx02),min(dy02));
mx = max(max(dx02),max(dy02));

%using rices_rule
% total_bins_sphere = 2*(length(dx02))^(1/3);%the bins contained in the spherical detector hits area
% total_bins_window = 4/pi*total_bins_sphere; %the bins contained within the detector window (approximately square)
% RES = floor(sqrt(total_bins_window)); %the total bins across the x or y dimension

%override
RES = 100;

Edges{1} = linspace(mi,mx,RES);
Edges{2} = Edges{1};

Ctr{1} = Edges{1} + (Edges{1}(2) - Edges{1}(1))/2;
Ctr{1}(end) = [];
Ctr{2} = Ctr{1};



N = histcounts2(dx02,dy02, Edges{1}, Edges{2});

%check what squaring signal does. 

%N = N.^2;

if sm == 1
    %N = smoothn(FDM,'robust');
    N = smoothn(N, 0.1);
end 

%% plot figure

%a -  y grid values
%b - x grid values

[b,a] = meshgrid(Ctr{1}, Ctr{2});

fig_dim_x = 6;
fig_dim_y = 5;

fig = figure('Units', 'centimeters');
fig.Position = [5,5,fig_dim_x, fig_dim_y];

h = surf(a,b, N);

%figure properties 
fsize = 9;
%colormap gray
%colormap jet
colormap inferno
%colormap parula
shading interp
view(2)
axis square
grid off
cb = colorbar;
set(get(cb,'title'), 'string', '{\rho} (ions/bin)', 'fontsize', fsize);
cb.FontSize = fsize;

%xlabel('Detector X(mm)','fontsize',fsize)
%ylabel('Detector Y(mm)','fontsize',fsize)
ax1 = gca;
%set(ax1,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on', 'xtick', -30:15:30,'ytick', -30:15:30)
set(ax1,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on')
axis off
xlim([min(dx) max(dx)])
ylim([min(dy) max(dy)])


%% change caxis

mid_c = median(N(N > 0));

std_c = std(N(N>0));

caxis([mid_c - 2*std_c, mid_c + 2*std_c]);

%% set paper properties
set(gcf, 'PaperPositionMode','auto');
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'color', 'white');
%axis off

end

