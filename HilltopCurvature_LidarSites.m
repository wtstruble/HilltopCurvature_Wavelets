%% Calculate hilltop curvature (wavelets & polynomial)

% This script walks you through how to extract catchments of interest and
% extract hilltop curvature, using both 2D polynomials and continuous
% wavelet transforms. While it is reasonably well automated, some sections
% do require manual input, so I recommend that when first using it that you
% go through section by section. The final sections focus on figure making,
% which you will possibly/likely find not relevant for your use. Adjust
% accordingly.

% The other main script in this directory, Synthetic.m, tests the wavelet
% and polynomial on synthetic hillslopes. 

% This script heavily utilizes TopoToolbox. Download @
% topotoolbox.wordpress.com

% Written by Will Struble, 2021, University of Arizona
% Published in Struble, W.T., Roering, J.J., 2021, Hilltop curvature as a
% proxy for erosion rate: Wavelets enable rapid computation and reveal
% systematic underestimation, Earth Surface Dynamics. 

%% Initial definitions
% An example DEM of an area for Hadsall Creek near the Siuslaw River in the Oregon
% Coast Range (Sius.tif) is included, as is a csv file with one set of
% coordinates (TestCoords.csv)

titlestring = 'Hadsall Creek'; % What is the name of this catchment?
DEM_path = 'Sius.tif';   %Input path to DEM 

erosion_file = 0;  % Do you have a csv file with erosion rates? If so, set to 1.
    % (Code was written with ==1 in mind. ==0 may have some issues).
coord_path_name = 'TestCoords.csv';    % If you have a file, enter the name here. Otherwise doesn't matter.

Z = GRIDobj(DEM_path);

%% 1) Curvature calculation (Make sure you are in the directory where you want
% curvature output files saved!)
% First we'll use wavelets, since they are quick! Do for both wavelet scale
% definitions. The input for the polynomial has to be odd, so we will only
% do odd values here for comparison. 
for i = 5:2:23  % For the Ricker wavelet, sigma >1, so i = 5 is the min (need odd value for poly)
    a_Las = i/(sqrt(2)*pi*Z.cellsize); % Lashermes et al. lambda definition
    a_Torr = (i*sqrt(5/2))/(2*pi*Z.cellsize); % Torrence and Compo definition
    
    C_Las = GRIDobj(Z);
    C_Torr = GRIDobj(Z);
    
    [C_Las.Z,~,~] = conv2_mexh_curv(Z.Z,a_Las,Z.cellsize);
    [C_Torr.Z,~,~] = conv2_mexh_curv(Z.Z,a_Torr,Z.cellsize);
    
    saveLas = ['C_Las_' num2str(i) 'm.tif'];
    saveTorr = ['C_Torr_' num2str(i) 'm.tif'];
    
    GRIDobj2geotiff(C_Las,saveLas);
    GRIDobj2geotiff(C_Torr,saveTorr);
    
    clear C_Las C_Torr a_Las a_Torr i saveLas saveTorr
    
end

% Now the slow part. Calculate curvatures with the polynomial for the same scales. 
for i = 5:2:23
    
    [~,~,C_poly] = poly2dgridv4(Z,i);
    savepoly = ['C_poly_' num2str(i) 'm.tif'];
    GRIDobj2geotiff(C_poly,savepoly);
    
    clear C_poly i savepoly 
    
end

% Now do the larger scales, but less often. We will also make this simpler
% and do it in one single for loop. 
for i = 25:4:45
    a_Las = i/(sqrt(2)*pi*Z.cellsize);
    a_Torr = (i*sqrt(5/2))/(2*pi*Z.cellsize);
    
    C_Las = GRIDobj(Z);
    C_Torr = GRIDobj(Z);
    
    [C_Las.Z,~,~] = conv2_mexh_curv(Z.Z,a_Las,Z.cellsize);
    [C_Torr.Z,~,~] = conv2_mexh_curv(Z.Z,a_Torr,Z.cellsize);
    
    saveLas = ['C_Las_' num2str(i) 'm.tif'];
    saveTorr = ['C_Torr_' num2str(i) 'm.tif'];
    
    GRIDobj2geotiff(C_Las,saveLas);
    GRIDobj2geotiff(C_Torr,saveTorr);
    
    clear C_Las C_Torr a_Las a_Torr saveLas saveTorr
    
    [~,~,C_poly] = poly2dgridv4(Z,i);
    savepoly = ['C_poly_' num2str(i) 'm.tif'];
    GRIDobj2geotiff(C_poly,savepoly);
    
    clear C_poly i savepoly 
    
end

% Now for the last ones at the biggest scales. Will take a while (underscored, in bold...).
for i = 51:10:141
    a_Las = i/(sqrt(2)*pi*Z.cellsize);
    a_Torr = (i*sqrt(5/2))/(2*pi*Z.cellsize);
    
    C_Las = GRIDobj(Z);
    C_Torr = GRIDobj(Z);
    
    [C_Las.Z,~,~] = conv2_mexh_curv(Z.Z,a_Las,Z.cellsize);
    [C_Torr.Z,~,~] = conv2_mexh_curv(Z.Z,a_Torr,Z.cellsize);
    
    saveLas = ['C_Las_' num2str(i) 'm.tif'];
    saveTorr = ['C_Torr_' num2str(i) 'm.tif'];
    
    GRIDobj2geotiff(C_Las,saveLas);
    GRIDobj2geotiff(C_Torr,saveTorr);
    
    clear C_Las C_Torr a_Las a_Torr saveLas saveTorr
    
    [~,~,C_poly] = poly2dgridv4(Z,i);
    savepoly = ['C_poly_' num2str(i) 'm.tif'];
    GRIDobj2geotiff(C_poly,savepoly);
    
    clear C_poly i savepoly 
    
end


% Finally, Gabet et al. use a smoothing scale of 14 m. We will do that for
% the wavelet transform for comparison (13 and 15 m is the closest we can
% get with the polynomial script we are using).
a_Las = 14/(sqrt(2)*pi*Z.cellsize);
a_Torr = (14*sqrt(5/2))/(2*pi*Z.cellsize);
    
C_Las = GRIDobj(Z);
C_Torr = GRIDobj(Z);
    
[C_Las.Z,~,~] = conv2_mexh_curv(Z.Z,a_Las,Z.cellsize);
[C_Torr.Z,~,~] = conv2_mexh_curv(Z.Z,a_Torr,Z.cellsize);
    
saveLas = ['C_Las_' num2str(i) 'm.tif'];
saveTorr = ['C_Torr_' num2str(i) 'm.tif'];
    
GRIDobj2geotiff(C_Las,saveLas);
GRIDobj2geotiff(C_Torr,saveTorr);
    
clear C_Las C_Torr a_Las a_Torr i saveLas saveTorr
    
%% 2) Import cosmogenic sampling locations
% (NOTE: THE CODE WILL ASSUME THAT THESE POINTS ARE CLOSE TO THE STREAM
% LOCATION WHERE THE SAMPLES WERE COLLECTED). IF YOU KNOW THEIR LOCATIONS
% ARE INACCURATE, YOU WILL NEED TO ADJUST THEM IN ARC FIRST (or do sectond
% option (i.e. erosion_file==0)).
% Assuming input is csv file and coordinates are still in lat/long. 
% Make sure csv columns are labeled Latitude and Longitude as well as
% ErosionRate and StandardDeviation.

if erosion_file==1
    
    T = readtable(coord_path_name);

% Assumes latitude is first column, longitude is second, Erosion rate is third:
    x = [T.Latitude];
    y = [T.Longitude];
    Sample_Erosion = [T.ErosionRate];
    Sample_Erosion_SD = [T.E_StandardDeviation];

else
    
% If you don't have sample sites, but want to manually select sites, then:
   % (Note this code wasn't originally written with this in mind, so there
   % may be some sections below (especially when figure making) that throw
   % some errors. Adjust as appropriate for your use). 
    hfig = figure;
    imageschs(Z);
    hold on
    % You may want to be informed by the stream network for selecting
    % points. If so, uncomment:  
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A = flowacc(FD).*(FD.cellsize^2);  % drainage area in m2
%     S = STREAMobj(FD,A>10000);
%     plot(S,'k');

    [x,y]=getpts(hfig);
    
end

%% 3) INCLUDES INTERACTIVE: Stream Networks
% Ok! So now that you have your curvatures for the whole rasters, we need
% to isolate the hilltops and clip the catchments where cosmogenic samples
% were collected. 

FD = FLOWobj(Z,'preprocess','carve');
S = STREAMobj(FD,flowacc(FD)>10000);   % The threshold drainage area
% shouldn't necessarily matter here, as the drainage basins are defined by
% the sample locations. 
% L = drainagebasins(FD);

if erosion_file==1
    hfig = figure; imageschs(Z); hold on; plot(S,'k');
    [x,y] = snap2stream(S,x,y,'inputislatlon',true);
else
    [x,y] = snap2stream(S,x,y); % If selected manually above.
end

xsave = x;
ysave = y;
Sample_Erosion_save = Sample_Erosion;
Sample_Erosion_SD_save = Sample_Erosion_SD;

%% 3a) Modify stream networks. NOTE!! This is where manual adjustment is necessary.
% Because I originally wrote this code to be run for a single basin at a
% time, you need to modify x and y and run them individually. I have this
% set to be done as automatically as I can. Rewriting the code to run
% through all sample sites is, as currently written, too memory intensive.

% Which sample in your csv (i.e. which row?) file are you running?
q = 1;    % Input row!
x = xsave(q);
y = ysave(q);
Sample_Erosion = Sample_Erosion_save(q);
Sample_Erosion_SD = Sample_Erosion_SD_save(q);

L = drainagebasins(FD,x,y);
Sc = modify(S,'clip',L);

% If your drainage network has other side tributaries that you don't need,
% manually clip out the main catchment(s) you want. Likely won't need
% this line. 
Sc = modify(Sc,'interactive','polyselect');

%% 4) Divides Mapping
% FYI, this may take some tinkering...

% Resample DEM to 5m; divides seem to plot better
Z_5m = resample(Z,5);
D = DIVIDEobj(FD,Sc,'outlets',false,'network',false);
D = divnet(D,FD);
D_s = sort(D);
D_o = divorder(D_s); %D_o = divdist(D_o);
D_f = removeshortdivides(D_o,FD,800);  % Remove divides shorter than 800 m
[x,y] = ind2coord(D_f,D_f.IX(:)); 

% Check: do the mapped divides look reasonable?
figure
imageschs(Z); hold on
plot(D_f,'color','k');

%% 5) Clean up divides (clip off edge effects) & convert to GRIDobj
bottomleftx_d = Z.refmat(3,1);
bottomrightx_d = bottomleftx_d + Z.cellsize*Z.size(2);
upperlefty_d = Z.refmat(3,2);
bottomlefty_d = upperlefty_d - Z.cellsize*Z.size(1);

Div = [x y]; 
Div(Div(:,1) < bottomleftx_d | Div(:,1) > bottomrightx_d) = NaN;
Div(Div(:,2) < bottomlefty_d | Div(:,2) > upperlefty_d) = NaN;

for i = 1:size(Div,1);
    if isnan(Div(i,1)) == 1 | isnan(Div(i,2)) == 1
        Div(i,:) = NaN;
    end
end

Divides = line2GRIDobj(Z,Div(:,1),Div(:,2));

clear bottomleftx_d bottomlefty_d bottomrightx_d upperlefty_d i

%% 6) Clip out low relief valleys (OPTIONAL! Will depend on landscape)
% (i.e. big rivers, alluvial valleys where there might be "divides" mapped))

H = resample(localtopography(resample(Z,5),50),Z); % Faster if you resample
H_log = (H<25); % Being fairly generous with this cutoff. Might (really depends) clip
% out some low-relief hilltops, but I think that is generally ok to make
% sure we don't have erroneous divides in the valleys. 
ix_log = find(H_log.Z==1);

Divides.Z(ix_log) = 0;

%% 7) Clip out main divide border
% In case there are any mobile divides (asymmetric) on the edges of our
% basin of interest, we will clip those out. (i.e. we will stick with internal divides)
% Could skip this step if not an issue in your landscape. 

[~,L_x,L_y] = GRIDobj2polygon(L);

L_GRID = line2GRIDobj(Z,L_x,L_y);

Divides.Z(L_GRID.Z==1) = 0;

clear L_GRID L_x L_y

%% 8) Final Clipping
% There can occasionally be erroneous divides that remain. This is
% particularly true if you are dealing with a landscape that has a lot of
% roads that may mess with flow routing. We want to clip these remaining
% divides out. 

% There is likely a better way to do this, but I have it written here as an
% iterative and manual (using the create mask) process, so you'll probably need to run these lines a few
% times. 

MASK = createmask(Divides); % Outline the divides that you DON'T want. 
% If you instead want to simply outline the divides you do want, change the
% next line to Divides.Z(MASK.Z==1) = 1;
Divides.Z(MASK.Z==1) = 0;

%% 8b) Run this after you're divide mask is correct.
% Save the mapped divide logical raster as a geotiff
GRIDobj2geotiff(Divides);

clear MASK

%% 9) INCLUDES INTERACTIVE: Clip out representative hilltop for finding scaling break (see Hurst et al., 2012; Struble and Roering, 2021)
MASK = createmask(Z,'usehillshade'); % You will be clipping on the blank hillshade, so make sure you know where the divide is. 
% Alternatively, you could use the logical raster, but I find seeing the
% hillshade useful for this step.
% Note what is "representative" is not included in this script. 

Divides_rep = Divides;
Divides_rep.Z(MASK.Z==0) = 0;

clear MASK

%% 10) Now clip curvatures by hilltops and find scaling break. 
fig1 = figure; hold on;

% First Lashermes
for i = 5:141
    file_read = strcat('C_Las_',string(i),'m.tif');
    
    if isfile(file_read) == 1
    Cht = GRIDobj(file_read);
    Cht_rep = Cht;
    
    Cht_rep.Z(Divides_rep.Z == 0) = NaN;
    
    fig1;
    subplot(3,1,1);
    hold on
    plot(i,nanmean(Cht_rep.Z,'all'),'ko','DisplayName','Lashermes et al.');
    subplot(3,1,2)
    hold on
    plot(i,nanstd(Cht_rep.Z,0,'all'),'ko','DisplayName','Lashermes et al.');
    subplot(3,1,3)
    hold on
    plot(i,iqr(Cht_rep.Z,'all'),'ko','DisplayName','Lashermes et al.');
      
    else
    end
end

% Now Torrence and Compo
for i = 5:141
    file_read = strcat('C_Torr_',string(i),'m.tif');
    
    if isfile(file_read) == 1
    Cht = GRIDobj(file_read);
    Cht_rep = Cht;
    
    Cht_rep.Z(Divides_rep.Z == 0) = NaN;
    
    fig1;
    subplot(3,1,1);
    hold on
    plot(i,nanmean(Cht_rep.Z,'all'),'rs','DisplayName','Torrence and Compo');
    subplot(3,1,2)
    hold on
    plot(i,nanstd(Cht_rep.Z,0,'all'),'rs','DisplayName','Torrence and Compo');
    subplot(3,1,3)
    hold on
    plot(i,iqr(Cht_rep.Z,'all'),'rs','DisplayName','Torrence and Compo');
    
    else
    end
end

for i = 5:141
    file_read = strcat('C_poly_',string(i),'m.tif');
    
    if isfile(file_read) == 1
    Cht = GRIDobj(file_read);
    Cht_rep = Cht;
    
    Cht_rep.Z(Divides_rep.Z == 0) = NaN;
    
    fig1;
    subplot(3,1,1);
    hold on
    plot(i,nanmean(Cht_rep.Z,'all'),'bd','DisplayName','Polynomial');
    subplot(3,1,2)
    hold on
    plot(i,nanstd(Cht_rep.Z,0,'all'),'bd','DisplayName','Polynomial');
    subplot(3,1,3)
    hold on
    plot(i,iqr(Cht_rep.Z,'all'),'bd','DisplayName','Polynomial');
  
    else
    end
end

% In case you want to see in log-scale...
subplot(3,1,1); set(gca,'xscale','log'); set(gca,'yscale','log'); title(titlestring); xlabel('Smoothing Scale (m)'); ylabel('C_{HT} (m^-^1)');
subplot(3,1,2); set(gca,'xscale','log'); set(gca,'yscale','log'); xlabel('Smoothing Scale (m)'); ylabel('St. Dev. C_{HT} (m^-^1)');
subplot(3,1,3); set(gca,'xscale','log'); set(gca,'yscale','log'); xlabel('Smoothing Scale (m)'); ylabel('Int. Qt. Rng. C_{HT} (m^-^1)');
legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial');

%% 11) Now compare the curvature values from each method to make sure they are consistent

% Will first try w/15 m scale
Cht_Las15 = GRIDobj('C_Las_15m.tif');
Cht_Torr15 = GRIDobj('C_Torr_15m.tif');
Cht_poly15 = GRIDobj('C_poly_15m.tif');

G = gradient8(Z); % Could use arcslope too. gradient8 is the more conservative approach,
                   % though, since it should likely clip out more.

% Define GRIDobjects
Cht_repLas15 = Cht_Las15; Cht_repTorr15 = Cht_Torr15; Cht_reppoly15 = Cht_poly15;

% First, on just the representative hilltops, with positive curvatures
% removed. Also remove spots where slope >0.4
Cht_repLas15.Z(Divides_rep.Z == 0) = NaN;    % Clip to indidvidual hilltop
Cht_repLas15.Z(Cht_repLas15.Z>=0) = NaN;     % Remove positive curvatures
Cht_repLas15.Z(G.Z>0.4) = NaN;                   % Remove high slopes

Cht_repTorr15.Z(Divides_rep.Z == 0) = NaN;    % Clip to indidvidual hilltop
Cht_repTorr15.Z(Cht_repTorr15.Z>=0) = NaN;     % Remove positive curvatures
Cht_repTorr15.Z(G.Z>0.4) = NaN;                   % Remove high slopes

Cht_reppoly15.Z(Divides_rep.Z == 0) = NaN;    % Clip to indidvidual hilltop
Cht_reppoly15.Z(Cht_reppoly15.Z>=0) = NaN;     % Remove positive curvatures
Cht_reppoly15.Z(G.Z>0.4) = NaN;                   % Remove high slopes

% Define PDFs
[pdf_repLas15,xpdf_repLas15] = ksdensity(Cht_repLas15.Z(isnan(Cht_repLas15.Z)==0),'Support',[-Inf 0],'BoundaryCorrection','reflection');
[pdf_repTorr15,xpdf_repTorr15] = ksdensity(Cht_repTorr15.Z(isnan(Cht_repTorr15.Z)==0),'Support',[-Inf 0],'BoundaryCorrection','reflection');
[pdf_reppoly15,xpdf_reppoly15] = ksdensity(Cht_reppoly15.Z(isnan(Cht_reppoly15.Z)==0),'Support',[-Inf 0],'BoundaryCorrection','reflection');

% Now all values in the catchment
Cht_Lasall15 = Cht_Las15;
Cht_Torrall15 = Cht_Torr15;
Cht_polyall15 = Cht_poly15;

Cht_Lasall15.Z(Divides.Z==0) = NaN;          % Remove non-divides
Cht_Torrall15.Z(Divides.Z==0) = NaN;
Cht_polyall15.Z(Divides.Z==0) = NaN;

Cht_Lasall15.Z(Cht_Lasall15.Z>=0) = NaN;        % Remove positive curvatures
Cht_Torrall15.Z(Cht_Torrall15.Z>=0) = NaN;
Cht_polyall15.Z(Cht_polyall15.Z>=0) = NaN;

Cht_Lasall15.Z(G.Z>0.4) = NaN;                   % Remove high slopes
Cht_Torrall15.Z(G.Z>0.4) = NaN;                   % Remove high slopes
Cht_polyall15.Z(G.Z>0.4) = NaN;                   % Remove high slopes

% PDFs
[pdf_allLas15,xpdf_allLas15] = ksdensity(Cht_Lasall15.Z(isnan(Cht_Lasall15.Z)==0),'Support',[-Inf 0],'BoundaryCorrection','reflection');
[pdf_allTorr15,xpdf_allTorr15] = ksdensity(Cht_Torrall15.Z(isnan(Cht_Torrall15.Z)==0),'Support',[-Inf 0],'BoundaryCorrection','reflection');
[pdf_allpoly15,xpdf_allpoly15] = ksdensity(Cht_polyall15.Z(isnan(Cht_polyall15.Z)==0),'Support',[-Inf 0],'BoundaryCorrection','reflection');


fig = figure;
figt = tiledlayout(1,1);
ax1 = axes(figt);
axis square

plot(ax1,xpdf_repLas15,pdf_repLas15,'k-'); hold on;
plot(ax1,xpdf_allLas15,pdf_allLas15,'k--');

plot(ax1,xpdf_repTorr15,pdf_repTorr15,'r-');
plot(ax1,xpdf_allTorr15,pdf_allTorr15,'r--');

plot(ax1,xpdf_reppoly15,pdf_reppoly15,'b-');
plot(ax1,xpdf_allpoly15,pdf_allpoly15,'b--');

sgtitle(titlestring);
xlabel('C_{HT} (m^-^1)')
ylabel('PDF')

legend('Lashermes et al., representative hilltop','Lashermes et al., all hilltops'...
    ,'Torrence and Compo, representative hilltop','Torrence and Compo, all hilltops',...
    'Polynomial, representative hilltop','Polynomial, all hilltops','Location','Northwest');

% And add erosion rate to the top x-axis
secylim = ylim;
secxlim = xlim;

ax2 = axes(figt);
axis square

Diffu = 0.003; % Transport Coefficient
E_range = linspace(0,1000*(-0.5)*Diffu*min(secxlim),10);
E_y = zeros(length(E_range),1);
plot(ax2,E_range,E_y,'k-','LineWidth',0.1);
set(ax2,'XDir','reverse');
xlabel(ax2,'Erosion Rate (mm yr^{-1})');
ylim(ax2,[0 secylim(2)]);

ax2.XAxisLocation = 'top';
ax2.Color = 'none';
ax1.Box = 'off'; clc

ax2.Box = 'off';
axis([ax1 ax2],'square');

%% 12) Calculate different erosion rates for each scale
Diffu = 0.003; % Transport Coefficient (NOTE!! This will vary between landscapes.
                % This is the value used for the Oregon Coast Range.)

E_15_Las = 1000*(-0.5)*Diffu*Cht_repLas15.Z;  % Representative Hilltops
E_15_Torr = 1000*(-0.5)*Diffu*Cht_repTorr15.Z;
E_15_poly = 1000*(-0.5)*Diffu*Cht_reppoly15.Z;

% Take the mean erosion rate
E_15_Las_mean = nanmean(E_15_Las,'all'); E_15_Las_StDev = nanstd(E_15_Las,0,'all');
E_15_Torr_mean = nanmean(E_15_Torr,'all'); E_15_Torr_StDev = nanstd(E_15_Torr,0,'all');
E_15_poly_mean = nanmean(E_15_poly,'all'); E_15_poly_StDev = nanstd(E_15_poly,0,'all');

E_15_Las_med= nanmedian(E_15_Las,'all');
E_15_Torr_med = nanmedian(E_15_Torr,'all');
E_15_poly_med = nanmedian(E_15_poly,'all');

%
E_15_Las_all = 1000*(-0.5)*Diffu*Cht_Lasall15.Z;  % All Hilltops
E_15_Torr_all = 1000*(-0.5)*Diffu*Cht_Torrall15.Z;
E_15_poly_all = 1000*(-0.5)*Diffu*Cht_polyall15.Z;

% Take the mean erosion rate. 
E_15_Las_allmean = nanmean(E_15_Las_all,'all'); E_15_Las_allStDev = nanstd(E_15_Las_all,0,'all');
E_15_Torr_allmean = nanmean(E_15_Torr_all,'all'); E_15_Torr_allStDev = nanstd(E_15_Torr_all,0,'all');
E_15_poly_allmean = nanmean(E_15_poly_all,'all'); E_15_poly_allStDev = nanstd(E_15_poly_all,0,'all');

% And median
E_15_Las_allmed = nanmedian(E_15_Las_all,'all'); 
E_15_Torr_allmed = nanmedian(E_15_Torr_all,'all');
E_15_poly_allmed = nanmedian(E_15_poly_all,'all');


figure(fig)
% a = xlim(ax1_had);
% xlimit = a(1)+0.1*(abs(a(1)));

xline(ax2,Sample_Erosion,'k-.');             % These values are all dependent on field site
xline(ax2,Sample_Erosion+Sample_Erosion_SD,'k:');
xline(ax2,Sample_Erosion-Sample_Erosion_SD,'k:');

xline(ax2,E_15_Las_allmean,'k--');
xline(ax2,E_15_Torr_allmean,'r--');
xline(ax2,E_15_poly_allmean,'b--');

xline(ax2,E_15_Las_mean,'k-');
xline(ax2,E_15_Torr_mean,'r-');
xline(ax2,E_15_poly_mean,'b-');

%% 13) Compare Methods as regression?
t = tiledlayout(3,1);

nexttile
plot(Cht_Lasall15.Z,Cht_polyall15.Z,'k.');
hold on
plot(Cht_repLas15.Z,Cht_reppoly15.Z,'r.');
xlabel('C_{HT} (m^{-1}), Lashermes et al.');
ylabel('C_{HT} (m^{-1}), Polynomial');
legend('All Hilltops','Representative Hilltop','Location','Northwest');
refline(1,0);
axis square; grid on;
plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], 'k-')
title(titlestring);

nexttile
plot(Cht_Torrall15.Z,Cht_polyall15.Z,'k.');
hold on
plot(Cht_repTorr15.Z,Cht_reppoly15.Z,'r.');
xlabel('C_{HT} (m^{-1}), Torrence and Compo');
ylabel('C_{HT} (m^{-1}), Polynomial');
refline(1,0);
axis square; grid on;
plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], 'k-')

nexttile
plot(Cht_Torrall15.Z,Cht_Lasall15.Z,'k.');
hold on
plot(Cht_repTorr15.Z,Cht_repLas15.Z,'r.');
xlabel('C_{HT} (m^{-1}), Torrence and Compo');
ylabel('C_{HT} (m^{-1}), Lashermes et al.');
refline(1,0);
axis square; grid on;
plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], 'k-')

t.Padding = 'none';
t.TileSpacing = 'none';
