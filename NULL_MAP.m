
function [N3,a,b] = NULL_MAP(dx, dy, nulls, sm)
%this function plots the change in electric field across the surface

%dx - inputted X co-ordinates
%dy - inputted Y co-ordinates
%nulls - the number of pulses between detection events
%sm - smooth SDM? (1 = yes, 0 = no)

%copyright (C) 2022
%author: A.J. Breen, The University of Sydney, 2022
%licence: BSD 2-Clause License - see LICENCE

fsize = 9;


%% calculating ratio for nulls

x1 = [];
y1 = [];
x2 = [];
y2 = [];


x1 = dx(nulls > 0 & nulls <= 10); %A small number of NULLS works well for zone lines.
y1 = dy(nulls > 0 & nulls <= 10);

%x1 = dx(nulls > 200); %greater than 200 sometimes works well for poles.
%y1 = dy(nulls > 200);

mi = min(min(x1),min(y1));
mx = max(max(x1),max(y1));

%using rices_rule
%RES = floor(2*length(x1)^(1/3));

%overide
RES = 100;
Edges{1} = linspace(mi,mx,RES);
Edges{2} = Edges{1};

Ctr{1} = Edges{1} + (Edges{1}(2) - Edges{1}(1))/2;
Ctr{1}(end) = [];
Ctr{2} = Ctr{1};

N1 = histcounts2(x1,y1, Edges{1}, Edges{2});

N3 = N1;


N3(isnan(N3))=0;
N3(isinf(N3))=0;

if sm == 1
    N3 = smoothn(N3, 0.1);
end 

[b,a] = meshgrid(Ctr{1}, Ctr{2});

%% plot figure

fig_dim_x = 6;
fig_dim_y = 5;

fig = figure('Units', 'centimeters');
fig.Position = [5,5,fig_dim_x, fig_dim_y];

h = surf(a,b, N3);

    colormap inferno
    shading interp
    view(2)
    axis square
    grid off
    cb = colorbar;
    set(get(cb,'title'), 'string', '{\rho} (ions/bin)', 'fontsize', 9);
    cb.FontSize = 9;
    %xlabel('Detector X(mm)','fontsize',fsize)
    %ylabel('Detector Y(mm)','fontsize',fsize)
    ax1 = gca;
    %set(ax1,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on', 'xtick', -30:15:30,'ytick', -30:15:30)
    set(ax1,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on')
    %axis off
    xlim([min(dx) max(dx)])
    ylim([min(dy) max(dy)])
    set(gcf, 'PaperPositionMode','auto');
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'color', 'white');
    axis off
    
%change caxis

mid_c = median(N3(N3 > 0));

std_c = std(N3(N3>0));

caxis([mid_c - 2*std_c, mid_c + 2*std_c]);
end     
     
 
 
 



