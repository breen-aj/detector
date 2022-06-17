
function [N_new,a,b] = EFIELD_MAP(dx, dy, m, sm, RANGES)
%this function plots the change in electric field across the surface

%dx - inputted X co-ordinates
%dy - inputted Y co-ordinates
%m - mass to charge ratio of ions in dataset
%sm - smooth SDM? (1 = yes, 0 = no)
%RANGES - the ranges corresponding to the different charge states of the
%species of interest

%copyright (C) 2022
%author: A.J. Breen, The University of Sydney, 2022
%licence: BSD 2-Clause License - see LICENCE

fsize = 12;


Range01 = RANGES(1,:);

Range02 = RANGES(2,:);


%% calculating ratio

numranges = size(Range01, 1);

x1 = [];
y1 = [];
x2 = [];
y2 = [];

% N1 = {};
% N2 = {};
% N3 = {};

for i = 1: numranges
x1 = [x1 dx(m > Range01(i, 1) & m < Range01(i, 2))];
y1 = [y1 dy(m > Range01(i, 1) & m < Range01(i, 2))];

x2 = [x2 dx(m > Range02(i, 1) & m < Range02(i, 2))];
y2 = [y2 dy(m > Range02(i, 1) & m < Range02(i, 2))];
    

%%calculating bin size
mi = min(min(dx),min(dy));
mx = max(max(dx),max(dy));

%using rices_rule
RES1 = floor(2*length(x1)^(1/3));
RES2 = floor(2*length(x2)^(1/3));

%take res as the average between two values
RES = floor(min(RES1,RES2)/2); %change to suit

Edges{1} = linspace(mi,mx,RES);
Edges{2} = Edges{1};

Ctr{1} = Edges{1} + (Edges{1}(2) - Edges{1}(1))/2;
Ctr{1}(end) = [];
Ctr{2} = Ctr{1};

end

N1 = histcounts2(x1,y1, Edges{1}, Edges{2}); 

N2 = histcounts2(x2,y2, Edges{1}, Edges{2});

%only consider bins where the count for each species is > 5
f1 = N1 > 5;
f2 = N2 > 5;

N1 = N1.*f1;
N2 = N2.*f2;

N3 = N1./N2;

N3(isnan(N3))=0;
N3(isinf(N3))=0;

if sm == 1
    N3 = smoothn(N3, 0.1);
end 

[b,a] = meshgrid(Ctr{1}, Ctr{2});

%% remove background hump

% normalise the plots and combine

%smooth the background

%N4 = smoothn(N3, 200);

Z = reshape(N3,[],1); %z co-ord
X = reshape(a,[],1); %x co-ord
Y = reshape(b,[],1); %y-co-ord

X = X(Z~=0);
Y = Y(Z~=0);
Z = Z(Z~=0);

%note: only fit to nice atoms in the centre (ideally removing background e.g. 10-20 mm out to avoid edge
%effects
rad_inner = 0;
rad_outer = 20;
X2 = X(X.^2 +Y.^2 < rad_outer.^2 & X.^2 +Y.^2 > rad_inner.^2 );
Y2 = Y(X.^2 +Y.^2 < rad_outer.^2 & X.^2 +Y.^2 > rad_inner.^2 );
Z2 = Z(X.^2 +Y.^2 < rad_outer.^2 & X.^2 +Y.^2 > rad_inner.^2 );


sf = fit([X2, Y2], Z2, 'poly22');
%fit is of the form Z = p00 + p10*x + p01*y + p20*x^2 + p11*x*y;
N_hump = sf.p00 + sf.p10*a + sf.p01*b + sf.p20*a.^2 + sf.p11*a.*b + sf.p02*b.^2; 

N_new = N3-N_hump;
N_new(N3 == 0) = 0;
N_new(N3 == 0) = min(min(N_new));
N_new = N_new-min(min(N_new));
%N_new = N_new.^2;

%% generate figure

fig_dim_x = 6; %change figure dimensions here
fig_dim_y = 5;

fig = figure('Units', 'centimeters');
fig.Position = [5,5,fig_dim_x, fig_dim_y];

h = surf(a,b, N_new);

%figure properties

    colormap inferno
    %colormap gray
    shading interp
    view(2)
    axis square
    grid off
    cb = colorbar;
    set(get(cb,'title'), 'string', 'Al^{+++}/Al^{++}', 'fontsize', 9);
    cb.FontSize = 9;
    ax1 = gca;
    set(ax1,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on')
    %axis off
    xlim([min(dx) max(dx)])
    ylim([min(dy) max(dy)])
    set(gcf, 'PaperPositionMode','auto');
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'color', 'white');
    axis off
  
mid_c = mean(N_new(N_new > 0));

std_c = std(N_new(N_new>0));

caxis([0, mid_c + 2*std_c]);
 
end     
     
 
 
 



