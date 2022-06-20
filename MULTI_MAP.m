
function [N1, a, b] = MULTI_MAP(dx, dy, Nat_pulse, sm)
%this function plots the change in electric field across the surface

%dx - inputted X co-ordinates
%dy - inputted Y co-ordinates
%Nat_pulse - the number of hits detected on pulse
%sm - smooth SDM? (1 = yes, 0 = no)


%copyright (C) 2022
%author: A.J. Breen, The University of Sydney, 2022
%licence: BSD 2-Clause License - see LICENSE

fsize = 9;

%% filter multiple hit ions

x1 = dx(Nat_pulse == 0 | Nat_pulse > 1); %note this filtering filters all multiple hit ions, a more sophisticated approach would be to specify the type of multiple hits to filter e.g. 2, 3.
y1 = dy(Nat_pulse == 0 | Nat_pulse > 1);


mi = min(min(x1),min(y1));
mx = max(max(y1),max(y1));

%using rices_rule
%RES = floor(2*length(x1)^(1/3));

%override
RES =100; 
Edges{1} = linspace(mi,mx,RES);
Edges{2} = Edges{1};

Ctr{1} = Edges{1} + (Edges{1}(2) - Edges{1}(1))/2;
Ctr{1}(end) = [];
Ctr{2} = Ctr{1};

%% calculating multiple hits

N1 = histcounts2(x1,y1, Edges{1}, Edges{2});

if sm == 1
    N1 = smoothn(N1, 0.1);
end 


[b,a] = meshgrid(Ctr{1}, Ctr{2});

%% plot figure

fig_dim_x = 6;%change the dimensions of the figure here
fig_dim_y = 5;

fig = figure('Units', 'centimeters');
fig.Position = [5,5,fig_dim_x, fig_dim_y];

h = surf(a,b, N1);

    colormap inferno
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
    %axis off
    xlim([min(dx) max(dx)])
    ylim([min(dy) max(dy)])

    set(gcf, 'PaperPositionMode','auto');
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'color', 'white');
    axis off
 
end     
     
 
 
 



