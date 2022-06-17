%SDM signal sweep of the reconstruction

close all

%count time to execute
tic


%open parallel pool
% % 
% if isempty(poolobj)
%         parpool
% end

poolobj = gcp;

%shut down parallel pool
%delete(poolobj)

%% user inputs
fsize = 9; %the font size
SDM_extent = 1;
SDM_bins_start = -SDM_extent;
SDM_bins_end = SDM_extent;

kf = 1.65; %default
ICF = 3.3;

% kf = 4.1; %Al dataset R18_53705
% ICF = 1.6;

% kf = 3.82; %Si dataset R18_60310
% ICF = 1.77;

% kf = 3.76; %Ti64 dataset R18_58152
% ICF = 1.7;

% 
% kf = 3.83; %Si dataset R04_09970
% ICF = 1.78;


% kf = 3.95; %Al dataset R18_53705
% ICF = 1.66;


% kf = 2.28; %Hansheng Cantor R18_57470
% ICF = 1.93;

% kf = 1.59; %Ni dataset R18_59205
% ICF = 3.84;

% kf = 3.3; %Ni dataset R18_59205
% ICF = 1.65;


% kf = 3.84; %In738 dataset R18_59377
% ICF = 1.59;

% kf = 3.93; %H282 dataset R04_23087
% ICF = 1.66;

avgDens = 91.4; %Ni
Fevap = 35; % V/nm from Miller textbook 

% avgDens = 60.66; %Al
% Fevap = 19; % V/nm from Miller textbook 

% avgDens = 85; %Fe
% Fevap = 33; % V/nm from Miller textbook 

%cantor alloy %Hansheng R18_57470
% avgDens = 84.7; %
% Fevap = 33; % V/nm from Miller textbook 

%si
% avgDens = 49.97; %
% Fevap = 33; % V/nm from Miller textbook 

%Ti
% avgDens = 57.27; %
% Fevap = 26; % V/nm from Miller textbook 

flightLength = 90;
detEff = 0.57;
%declare percentage array

global perc %make global
perc  = zeros(1,10); 
%% read the epos
filename = 'R18_59205-v01.epos';
[dx,dy,~,~,~,m,~,vdc,vp,~,~] = readepos(filename);

SaveFileName = [filename 'signal_map_all02.png'];

%% region of interest

%radius of ROI
R = 2;

%ion number limits (z slice)
s1 = 3e6;
s2 = 5e6;

%%find extents of the data

dx02 = dx(s1:s2);
dy02 = dy(s1:s2);
m02 = m(s1:s2);
vdc02 = vdc(s1:s2);
vp02 = vp(s1:s2);

%detector area crop
rad = 30;
keep = dx02.^2 + dy02.^2 < rad.^2;

dx02 = dx02(keep);
dy02 = dy02(keep);
m02 = m02(keep);
vdc02 = vdc02(keep);
vp02 = vp02(keep);

% x,y extents of slice

dxmin = min(dx02)-0.01;
dxmax = max(dx02)+0.01;

dymin = min(dy02)-0.01;
dymax = max(dy02)+0.01;

%sample spacing

s_space = 1; %mm %note use the interp2 function to change the spacing/binning to suit other maps. 

xbins = ceil(norm(dxmax-dxmin)/s_space)+1;
ybins = ceil(norm(dymax-dymin)/s_space)+1;

x_mid = (dxmax+dxmin)/2;
y_mid = (dymax+dymin)/2;

x_start = floor(x_mid - s_space*xbins/2);
y_start = floor(y_mid - s_space*ybins/2);

[X,Y] = meshgrid(x_start:s_space: x_start + xbins*s_space, y_start:s_space: y_start + ybins*s_space); %edges of voxels

X_cent = X +s_space/2; %centres of voxels
X_cent(:, end) = [];
X_cent(end, :) = [];

Y_cent = Y +s_space/2; %centres of voxels
Y_cent(:, end) = [];
Y_cent(end, :) = [];


%SDM extents
%smallest bin size of IVAS SDM = 0.008 nm
edges = SDM_bins_start:0.008:SDM_bins_end;
bins = length(edges)-1;

centres = edges+(edges(2)-edges(1))/2;
centres(end) = [];

%total number of bins
tot = (xbins)*(ybins);

count = 0;
Int_val = zeros(ybins, xbins);

nLoops = xbins;

updateWaitbar = waitbarParfor(nLoops, "Calculation in progress...");

%calculate if voxel is on edge
[N, x_edge_test, y_edge_test] = histcounts2(dx02, dy02, X(1, :), Y(:, 1));
N = N';
parfor i = 1:xbins
    updateWaitbar();
    for j = 1:ybins
        
%         if j == 2 & i == 5
%             disp('here');
%         end
        
        if N(j, i) < 10 %set this value based on the
            Int_val(j, i) = 0;
            continue
        end
        h = x_start+s_space*(i-1);
        k = y_start+s_space*(j-1);
        
        %override
        %P1
%         h = 5.024;
%         k = 9.559;
%         
        %noise
%         h = -2.004;
%         k = -4.952;
        
        l1 = (dx02-h).^2 + (dy02-k).^2 < R.^2;
        dx03 = dx02(l1);
%         if i == 24 && j ==24
%             disp('stop');
%         end
        
        if isempty(dx03)
            Int_val(j,i) = 0;
            count = count +1;
            %disp([num2str(count) ' SDMs calculated of ' num2str(tot)])
            continue
        end
        dy03 = dy02(l1);
        vdc03 = vdc02(l1);
        vp03 = vp02(l1);
        m03 = m02(l1);
        
        
        %[x, y, z, x2, y2, z2] = atomProbeRecon_SDM(dx03, dy03, h, k, vdc03+vp03, kf, ICF, SDM_extent);
        %[x, y, z, m03, x2, y2, z2, m04] = atomProbeRecon_SDM(dx03, dy03, m03,h,k, vdc03+vp03, kf, ICF, 90, 0.57, avgDens, Fevap, SDM_extent);
        [x, y, z, R_start, R_end] = atomProbeRecon05(dx03, dy03,h,k, vdc03+vp03, kf, ICF, avgDens, Fevap, flightLength, detEff);
        %% compute SDM for species of interest
        
        filt_element = true;
        
        if filt_element == true
        
        
        %RANGES = [8.847, 9.146; 13.372, 13.898]; %Al
        RANGES = [13.375, 13.818]; %Al
        %RANGES = [13.347, 13.898; 26.855, 27.357]; %2+, 1+
       
        
        x_filt = [];
        y_filt = [];
        z_filt = [];
        m_filt = [];
        for i2 = 1:size(RANGES, 1)
           x_filt = [x_filt x(m03 > RANGES(i2,1) & m03 < RANGES(i2,2))];
           y_filt = [y_filt y(m03 > RANGES(i2,1) & m03 < RANGES(i2,2))];
           z_filt = [z_filt z(m03 > RANGES(i2,1) & m03 < RANGES(i2,2))];
        end
       
        else
            z_filt = z;
        end

%         if count == 73
%             disp('stop');
%         end
        sdm = zeros(1, bins);
        
%     if isempty(x2) == 0    
        for j2 = 1:length(z_filt)
            dz = z_filt - z_filt(j2);
            k = dz > SDM_bins_start & dz < SDM_bins_end;
            sdm = sdm + histcounts(dz(k), edges);
        end
%     else
%         for j2 = 1:length(z) %this is if thickness of slice is too small for inmputed SDM extents
%             dz = z - z(j2);
%             k = dz > SDM_bins_start & dz < SDM_bins_end;
%             sdm2 = sdm2 + histcounts(dz(k), edges);
%         end
%     end
 
   %remove atoms counting themselves
   %sdm2((length(edges)+1)/2) = sdm2((length(edges)+1)/2) - length(z2);
%   Int_val(j, i) = var(sdm2); %index based purely on variace

%make symmetric since in this case not all atoms are counting themselves twice
%sdm2 = sdm2 +flip(sdm2);

%% plot the SDM

% centres = edges+(edges(2)-edges(1))/2;
% centres(end) = [];

% plot(centres, sdm, 'linewidth', 2, 'color', 'b');
% 
% 
% set(gca,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on');
% 
% xlim([SDM_bins_start SDM_bins_end]);
% 
% xlabel('z'' (nm)', 'fontsize', fsize);
% ylabel('counts', 'fontsize', fsize);
% 
% set(gcf, 'PaperPositionMode','auto');
% set(gcf, 'InvertHardCopy', 'off');
% set(gcf, 'color', 'white');
% 
% legend('SDM');


%% normalised var measure
sdm_sum = sum(sdm);
%if sdm_sum > 100 %only use cells where enough atoms are present to get a reasonable SDM i.e. 100 would mean 10 atoms in the voxel atoms^2 = sdm_sum. 
   %Int_val(j, i) = var(sdm2)/(sum(sdm2)); %normalised to number of counts in SDM
   %calculate metric based on DF fit. 
   P = polyfit(centres, sdm, 2); %polynomial fit of degree 2
   F_fit = P(1)*centres.^2 + P(2)*centres + P(3);
   %Int_val(j, i) = 1/(length(z).^2)*trapz((sdm2 - F_fit).^2); 
   Int_val(j,i) = trapz((sdm - F_fit).^2)*1/sdm_sum;
   %Int_val(j,i) = 1/(length(z)*length(z2))*trapz((sdm2 - F_fit).^2);
   %Int_val(j,i) = 1/(sum(sdm2))*trapz((sdm2 - F_fit).^2);
   
% else
%    Int_val(j, i) = 0;
% end
   
   
   
   %Int_val(j, i) = var(sdm2)*length(x);  %be careful, remember first position refers to rows, which represents y and second position refers to columns, which represents x. i.e. (y, x) This is different to the usual co-ordinate convention (x,y)
   %weight the variance by the number of atoms in the ROI
   %count = count +1;
   
   % calculate % complete
   %progress_percent(count,tot)
   
   %disp([num2str(count) ' SDMs calculated of ' num2str(tot)])
   
    end
end

%% plot figure

Int_val(isnan(Int_val))=0;
Int_val(isinf(Int_val))=0;

fig_dim_x = 5;
fig_dim_y = 4;

fig = figure('Units', 'centimeters');
fig.Position = [5,5,fig_dim_x, fig_dim_y];

surf(X_cent, Y_cent, Int_val);

colormap inferno
%colormap gray
shading interp
view(2)
axis square
grid off
b = colorbar;
set(get(b,'title'), 'string', 'E', 'fontsize', 9);
b.FontSize = 9;

%change caxis of colourmap

mid_c = median(Int_val(Int_val > 0));

std_c = std(Int_val(Int_val>0));

%caxis([mid_c - 2*std_c, mid_c + 2*std_c]); 
caxis([0, mid_c + 2*std_c]);
%xlabel(' det X (mm)','fontsize',fsize)
%ylabel(' det Y (mm)','fontsize',fsize)
ax1 = gca;
set(ax1,'XColor','k','YColor','k','fontsize',fsize,'linewidth', 2, 'box','on')
axis off
xlim([X_cent(1,1) X_cent(1, end)]);
ylim([Y_cent(1,1) Y_cent(end, 1)])
set(gcf, 'PaperPositionMode','auto');
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'color', 'white');     
export_fig(SaveFileName, '-png', '-r500');

toc       
        

