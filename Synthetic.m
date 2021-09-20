%% Defining a synthetic hillslope with constant curvature to compare curvature calculation methods

% This script constructs synthetic hillslopes and calculates hilltop
% curvature using 2D polynomials and continuous wavelet transforms. It then
% compares the outputs and determines whether hilltop curvature is
% underestimated as a function of erosion rate and/or smoothing scale. 

% Written by Will Struble, 2021, University of Arizona
% Published in Struble, W.T., Roering, J.J., 2021, Hilltop curvature as a
% proxy for erosion rate: Wavelets enable rapid computation and reveal
% systematic underestimation, Earth Surface Dynamics. 


%% Defining a 1-D nonlinear diffusive hillslope profile translated in the y-direction

[xspace_t,yspace_t] = meshgrid(-100:100,-100:100,1);  % Make 201x201 matrix
                    % (201 so that we have a middle pixel exactly at the hilltop

%% Apply nonlinear diffusive theoretical profile from Roering et al.(2007)

K = 0.003; % This is the diffusivity. Using Oregon Coast Range value here. 
Lh = 101;  % Hillslope length.
           % 201x201 grid (there would be no middle row if even dimensions)
Sc125 = 1.25;  % Critical slope angle
Sc = Sc125; % Could use different Sc

Estar1 = 1;  % Setting up E* values to be considered
Estar10 = 10;
Estar30 = 30;
Estar100 = 100;

Cht1 = (Estar1*Sc)/(2*Lh);
Cht10 = (Estar10*Sc)/(2*Lh);
Cht30 = (Estar30*Sc)/(2*Lh);
Cht100 = (Estar100*Sc)/(2*Lh);

eros1 = K*Cht1*0.5;
eros10 = K*Cht10*0.5;
eros30 = K*Cht30*0.5;
eros100 = K*Cht100*0.5;

% IMPORTANT: To make sure we are using the same random numbers through
% all tests, we will NOT rerun the following lines when we test the other
% parameters (if we do redefine R, Rp, or Rr, we will need to rerun all the
% parameter tests, to make sure we are consistent). 

% White Noise
R = randnd(0,201);  % randnd produces normally distributed values with 1sigma between -1 and 1 (see Konowalczyk (2021) on FileExchange)
% Pink noise
Rp = randnd(-1,201); 
% Red noise
Rr = randnd(-2,201); 

% Apply the wavelet and polynomial to scales from 5m to 35 m (odd numbers).
% First define empty cells. Note that this will be a 3-dimensional cell
% array (m x j x iter = white/pink/red noise x smoothing scale x iteration)
C_Las_white = {}; C_Torr_white = {}; C_poly_white = {};
C_Las_pink = {}; C_Torr_pink = {}; C_poly_pink = {}; 
C_Las_red = {}; C_Torr_red = {}; C_poly_red = {};


Estar_Las_white = {}; Estar_Torr_white = {}; Estar_poly_white = {}; 
Estar_Las_pink = {}; Estar_Torr_pink = {}; Estar_poly_pink = {};
Estar_Las_red = {}; Estar_Torr_red = {}; Estar_poly_red = {};

%% First iteration: Very low, low, moderate, high E*, change beta; sigma = 0.5% Lh
iter = 1;  % Which iteration (i.e. noise amplitude) is this?

z_Roer1 = ((K*Sc125^2)/(2*2*eros1))*(log(0.5*(sqrt(1+((2*2*eros1*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros1*xspace_t)/(K*Sc125)).^2)+1);
z_Roer10 = ((K*Sc125^2)/(2*2*eros10))*(log(0.5*(sqrt(1+((2*2*eros10*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros10*xspace_t)/(K*Sc125)).^2)+1);
z_Roer30 = ((K*Sc125^2)/(2*2*eros30))*(log(0.5*(sqrt(1+((2*2*eros30*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros30*xspace_t)/(K*Sc125)).^2)+1);
z_Roer100 = ((K*Sc125^2)/(2*2*eros100))*(log(0.5*(sqrt(1+((2*2*eros100*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros100*xspace_t)/(K*Sc125)).^2)+1);

% Noise Amplitude
Amp = 0.005*Lh;   % Noise Amplitude is 0.5% of hillslope length 

% White Noise
z_Roer1_R = z_Roer1 + Amp*R;  % Noise matrix
z_Roer10_R = z_Roer10 + Amp*R;
z_Roer30_R = z_Roer30 + Amp*R;
z_Roer100_R = z_Roer100 + Amp*R;


z_Roer1_R = z_Roer1_R + abs(nanmin(z_Roer1_R,[],'all'));    % Reset baselevel as 0
z_Roer10_R = z_Roer10_R + abs(nanmin(z_Roer10_R,[],'all'));    % Reset baselevel as 0
z_Roer30_R = z_Roer30_R + abs(nanmin(z_Roer30_R,[],'all'));       % Reset baselevel as 0
z_Roer100_R = z_Roer100_R + abs(nanmin(z_Roer100_R,[],'all'));       % Reset baselevel as 0


% Pink Noise
z_Roer1_Rp = z_Roer1 + Amp*Rp;
z_Roer10_Rp = z_Roer10 + Amp*Rp;
z_Roer30_Rp = z_Roer30 + Amp*Rp;
z_Roer100_Rp = z_Roer100 + Amp*Rp;


z_Roer1_Rp = z_Roer1_Rp + abs(nanmin(z_Roer1_Rp,[],'all'));    % Reset baselevel as 0
z_Roer10_Rp = z_Roer10_Rp + abs(nanmin(z_Roer10_Rp,[],'all'));    % Reset baselevel as 0
z_Roer30_Rp = z_Roer30_Rp + abs(nanmin(z_Roer30_Rp,[],'all'));       % Reset baselevel as 0
z_Roer100_Rp = z_Roer100_Rp + abs(nanmin(z_Roer100_Rp,[],'all'));       % Reset baselevel as 0


% Red Noise
z_Roer1_Rr = z_Roer1 + Amp*Rr;
z_Roer10_Rr = z_Roer10 + Amp*Rr;
z_Roer30_Rr = z_Roer30 + Amp*Rr;
z_Roer100_Rr = z_Roer100 + Amp*Rr;

z_Roer1_Rr = z_Roer1_Rr + abs(nanmin(z_Roer1_Rr,[],'all'));    % Reset baselevel as 0
z_Roer10_Rr = z_Roer10_Rr + abs(nanmin(z_Roer10_Rr,[],'all'));    % Reset baselevel as 0
z_Roer30_Rr = z_Roer30_Rr + abs(nanmin(z_Roer30_Rr,[],'all'));       % Reset baselevel as 0
z_Roer100_Rr = z_Roer100_Rr + abs(nanmin(z_Roer100_Rr,[],'all'));       % Reset baselevel as 0


j = 0;
for i=5:2:35
    j = j+1;
    a_Las_syn = i/(sqrt(2)*pi*1);
    a_Torr_syn = (i*sqrt(5/2))/(2*pi*1);
    
    % white noise
    [C_Las_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Las_syn,1);
    [C_Torr_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR1] = poly2dgridv4_synth(z_Roer1_R,i);
    
    m = 1; % First row is very low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR1;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR1;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR1;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR1.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR1.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR1.*Lh)./Sc;
    
    [C_Las_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Las_syn,1);
    [C_Torr_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR10] = poly2dgridv4_synth(z_Roer10_R,i);
    
    m = 2; % Second row is low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR10;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR10;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR10;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR10.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR10.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Las_syn,1);
    [C_Torr_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR30] = poly2dgridv4_synth(z_Roer30_R,i);
    
    m = 3; % Third row is sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR30;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR30;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR30;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR30.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR30.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR30.*Lh)./Sc;
   
    
    [C_Las_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Las_syn,1);
    [C_Torr_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR100] = poly2dgridv4_synth(z_Roer100_R,i);
    
    m = 4; % Third row is very sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR100;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR100;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR100;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR100.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR100.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR100.*Lh)./Sc;
    
    % pink noise
    [C_Las_syn_RoerRp1] = conv2_mexh_curv(z_Roer1_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp1] = conv2_mexh_curv(z_Roer1_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp1] = poly2dgridv4_synth(z_Roer1_Rp,i);
    
    m = 1; % First row is very low curvature hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp1;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp1;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp1;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp1.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp1.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp1.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp10] = conv2_mexh_curv(z_Roer10_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp10] = conv2_mexh_curv(z_Roer10_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp10] = poly2dgridv4_synth(z_Roer10_Rp,i);
    
    m = 2; % Second row is gentle hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp10;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp10;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp10;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp10.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp10.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp30] = conv2_mexh_curv(z_Roer30_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp30] = conv2_mexh_curv(z_Roer30_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp30] = poly2dgridv4_synth(z_Roer30_Rp,i);
    
    m = 3; % Third row is sharp hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp30;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp30;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp30;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp30.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp30.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp30.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp100] = conv2_mexh_curv(z_Roer100_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp100] = conv2_mexh_curv(z_Roer100_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp100] = poly2dgridv4_synth(z_Roer100_Rp,i);
    
    m = 4; % Third row is sharp hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp100;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp100;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp100;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp100.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp100.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp100.*Lh)./Sc;
    
    
    % red noise
    [C_Las_syn_RoerRr1] = conv2_mexh_curv(z_Roer1_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr1] = conv2_mexh_curv(z_Roer1_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr1] = poly2dgridv4_synth(z_Roer1_Rr,i);
    
    m = 1; % First row is moderate curvature hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr1;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr1;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr1;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr1.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr1.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr1.*Lh)./Sc;
    
    [C_Las_syn_RoerRr10] = conv2_mexh_curv(z_Roer10_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr10] = conv2_mexh_curv(z_Roer10_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr10] = poly2dgridv4_synth(z_Roer10_Rr,i);
    
    m = 2; % Second row is gentle hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr10;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr10;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr10;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr10.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr10.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRr30] = conv2_mexh_curv(z_Roer30_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr30] = conv2_mexh_curv(z_Roer30_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr30] = poly2dgridv4_synth(z_Roer30_Rr,i); 
    
    m = 3; % Third row is sharp hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr30;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr30;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr30;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr30.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr30.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr30.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRr100] = conv2_mexh_curv(z_Roer100_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr100] = conv2_mexh_curv(z_Roer100_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr100] = poly2dgridv4_synth(z_Roer100_Rr,i); 
    
    m = 4; % Third row is sharp hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr100;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr100;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr100;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr100.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr100.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr100.*Lh)./Sc;
    
end

%% First iteration: Plot Results, High, med, low Cht, change beta
figure
hold on
j = 0; iter = 1;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the moderate curvature 
    subplot(4,3,1)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,2)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,3)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,5)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,6)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,8)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,9)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

    
    m = 4; % Now the very sharp hilltops 
    subplot(4,3,10)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,11)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,12)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

end

subplot(4,3,1); yline((-1)*Cht1,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline((-1)*Cht1,'k--'); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline((-1)*Cht1,'k--'); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline((-1)*Cht10,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,5); yline((-1)*Cht10,'k--'); xlim([0 40])
subplot(4,3,6); yline((-1)*Cht10,'k--'); xlim([0 40])
subplot(4,3,7); yline((-1)*Cht30,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,8); yline((-1)*Cht30,'k--'); xlim([0 40])
subplot(4,3,9); yline((-1)*Cht30,'k--'); xlim([0 40])
subplot(4,3,10); yline((-1)*Cht100,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline((-1)*Cht100,'k--'); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline((-1)*Cht100,'k--'); xlabel('Smoothing scale(m)','FontSize',13); xlim([0 40])



legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('C_{HT}, \sigma = 0.5% L_H','FontSize',20);

%% First iteration: Plot ratio of calculated to known E*
figure
hold on
j = 0; iter = 1;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the very low curvatures 
    subplot(4,3,1)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    subplot(4,3,2)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    subplot(4,3,3)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    subplot(4,3,5)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    subplot(4,3,6)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    subplot(4,3,8)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    subplot(4,3,9)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    
    m = 4; % Now the sharp hilltops 
    subplot(4,3,10)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');
    
    subplot(4,3,11)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');
    
    subplot(4,3,12)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');

end

subplot(4,3,1); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline(1); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline(1); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlim([0 40])
subplot(4,3,5); yline(1); xlim([0 40])
subplot(4,3,6); yline(1); xlim([0 40])
subplot(4,3,7); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlim([0 40])
subplot(4,3,8); yline(1); xlim([0 40])
subplot(4,3,9); yline(1); xlim([0 40])
subplot(4,3,10); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline(1); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline(1); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])

legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('C_{HT} Ratio, \sigma = 0.5% L_H','FontSize',20);

%% First iteration: Plot standard deviation of E*
figure
hold on
j = 0; iter = 1;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the very low curvatures 
    subplot(4,3,1)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,2)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,3)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,5)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,6)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,8)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,9)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 4; % Now the sharp hilltops 
    subplot(4,3,10)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,11)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,12)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

end

subplot(4,3,1); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline(0); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline(0); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,5); yline(0); xlim([0 40])
subplot(4,3,6); yline(0); xlim([0 40])
subplot(4,3,7); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,8); yline(0); xlim([0 40])
subplot(4,3,9); yline(0); xlim([0 40])
subplot(4,3,10); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline(0); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline(0); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])

legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('Standard Deviation C_{HT}, \sigma = 0.5% L_H','FontSize',20);

%% Second Iteration: A = 0.1%

iter = 2;  % Which iteration is this?

% Noise Amplitude
Amp = 0.001*Lh;   % Noise Amplitude is 1% of hillslope length 

% White Noise
z_Roer1_R = z_Roer1 + Amp*R;  % Noise matrix
z_Roer10_R = z_Roer10 + Amp*R;
z_Roer30_R = z_Roer30 + Amp*R;
z_Roer100_R = z_Roer100 + Amp*R;


z_Roer1_R = z_Roer1_R + abs(nanmin(z_Roer1_R,[],'all'));    % Reset baselevel as 0
z_Roer10_R = z_Roer10_R + abs(nanmin(z_Roer10_R,[],'all'));    % Reset baselevel as 0
z_Roer30_R = z_Roer30_R + abs(nanmin(z_Roer30_R,[],'all'));       % Reset baselevel as 0
z_Roer100_R = z_Roer100_R + abs(nanmin(z_Roer100_R,[],'all'));       % Reset baselevel as 0


% Pink Noise
z_Roer1_Rp = z_Roer1 + Amp*Rp;
z_Roer10_Rp = z_Roer10 + Amp*Rp;
z_Roer30_Rp = z_Roer30 + Amp*Rp;
z_Roer100_Rp = z_Roer100 + Amp*Rp;


z_Roer1_Rp = z_Roer1_Rp + abs(nanmin(z_Roer1_Rp,[],'all'));    % Reset baselevel as 0
z_Roer10_Rp = z_Roer10_Rp + abs(nanmin(z_Roer10_Rp,[],'all'));    % Reset baselevel as 0
z_Roer30_Rp = z_Roer30_Rp + abs(nanmin(z_Roer30_Rp,[],'all'));       % Reset baselevel as 0
z_Roer100_Rp = z_Roer100_Rp + abs(nanmin(z_Roer100_Rp,[],'all'));       % Reset baselevel as 0


% Red Noise
z_Roer1_Rr = z_Roer1 + Amp*Rr;
z_Roer10_Rr = z_Roer10 + Amp*Rr;
z_Roer30_Rr = z_Roer30 + Amp*Rr;
z_Roer100_Rr = z_Roer100 + Amp*Rr;

z_Roer1_Rr = z_Roer1_Rr + abs(nanmin(z_Roer1_Rr,[],'all'));    % Reset baselevel as 0
z_Roer10_Rr = z_Roer10_Rr + abs(nanmin(z_Roer10_Rr,[],'all'));    % Reset baselevel as 0
z_Roer30_Rr = z_Roer30_Rr + abs(nanmin(z_Roer30_Rr,[],'all'));       % Reset baselevel as 0
z_Roer100_Rr = z_Roer100_Rr + abs(nanmin(z_Roer100_Rr,[],'all'));       % Reset baselevel as 0


j = 0;
for i=5:2:35
    j = j+1;
    a_Las_syn = i/(sqrt(2)*pi*1);
    a_Torr_syn = (i*sqrt(5/2))/(2*pi*1);
    
    % white noise
    [C_Las_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Las_syn,1);
    [C_Torr_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR1] = poly2dgridv4_synth(z_Roer1_R,i);
    
    m = 1; % First row is very low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR1;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR1;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR1;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR1.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR1.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR1.*Lh)./Sc;
    
    
    [C_Las_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Las_syn,1);
    [C_Torr_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR10] = poly2dgridv4_synth(z_Roer10_R,i);
    
    m = 2; % Second row is low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR10;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR10;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR10;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR10.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR10.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Las_syn,1);
    [C_Torr_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR30] = poly2dgridv4_synth(z_Roer30_R,i);
    
    m = 3; % Third row is sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR30;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR30;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR30;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR30.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR30.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR30.*Lh)./Sc;
   
    
    [C_Las_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Las_syn,1);
    [C_Torr_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR100] = poly2dgridv4_synth(z_Roer100_R,i);
    
    m = 4; % Third row is very sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR100;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR100;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR100;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR100.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR100.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR100.*Lh)./Sc;
    
    
    % pink noise
    [C_Las_syn_RoerRp1] = conv2_mexh_curv(z_Roer1_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp1] = conv2_mexh_curv(z_Roer1_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp1] = poly2dgridv4_synth(z_Roer1_Rp,i);
    
    m = 1; % First row is very low curvature hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp1;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp1;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp1;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp1.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp1.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp1.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp10] = conv2_mexh_curv(z_Roer10_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp10] = conv2_mexh_curv(z_Roer10_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp10] = poly2dgridv4_synth(z_Roer10_Rp,i);
    
    m = 2; % Second row is gentle hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp10;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp10;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp10;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp10.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp10.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp30] = conv2_mexh_curv(z_Roer30_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp30] = conv2_mexh_curv(z_Roer30_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp30] = poly2dgridv4_synth(z_Roer30_Rp,i);
    
    m = 3; % Third row is sharp hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp30;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp30;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp30;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp30.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp30.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp30.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp100] = conv2_mexh_curv(z_Roer100_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp100] = conv2_mexh_curv(z_Roer100_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp100] = poly2dgridv4_synth(z_Roer100_Rp,i);
    
    m = 4; % Fourth row is very sharp hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp100;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp100;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp100;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp100.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp100.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp100.*Lh)./Sc;
    
    
    % red noise
    [C_Las_syn_RoerRr1] = conv2_mexh_curv(z_Roer1_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr1] = conv2_mexh_curv(z_Roer1_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr1] = poly2dgridv4_synth(z_Roer1_Rr,i);
    
    m = 1; % First row is moderate curvature hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr1;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr1;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr1;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr1.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr1.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr1.*Lh)./Sc;
    
    [C_Las_syn_RoerRr10] = conv2_mexh_curv(z_Roer10_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr10] = conv2_mexh_curv(z_Roer10_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr10] = poly2dgridv4_synth(z_Roer10_Rr,i);
    
    m = 2; % Second row is gentle hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr10;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr10;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr10;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr10.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr10.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRr30] = conv2_mexh_curv(z_Roer30_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr30] = conv2_mexh_curv(z_Roer30_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr30] = poly2dgridv4_synth(z_Roer30_Rr,i); 
    
    m = 3; % Third row is sharp hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr30;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr30;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr30;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr30.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr30.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr30.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRr100] = conv2_mexh_curv(z_Roer100_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr100] = conv2_mexh_curv(z_Roer100_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr100] = poly2dgridv4_synth(z_Roer100_Rr,i); 
    
    m = 4; % Third row is sharp hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr100;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr100;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr100;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr100.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr100.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr100.*Lh)./Sc;
    
end

%% Second iteration: Plot Results, High, med, low Cht, change beta
figure
hold on
j = 0; iter = 2;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the moderate curvature 
    subplot(4,3,1)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,2)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,3)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,5)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,6)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,8)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,9)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

    
    m = 4; % Now the very sharp hilltops 
    subplot(4,3,10)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,11)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,12)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

end

subplot(4,3,1); yline((-1)*Cht1,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline((-1)*Cht1,'k--'); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline((-1)*Cht1,'k--'); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline((-1)*Cht10,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,5); yline((-1)*Cht10,'k--'); xlim([0 40])
subplot(4,3,6); yline((-1)*Cht10,'k--'); xlim([0 40])
subplot(4,3,7); yline((-1)*Cht30,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,8); yline((-1)*Cht30,'k--'); xlim([0 40])
subplot(4,3,9); yline((-1)*Cht30,'k--'); xlim([0 40])
subplot(4,3,10); yline((-1)*Cht100,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline((-1)*Cht100,'k--'); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline((-1)*Cht100,'k--'); xlabel('Smoothing scale(m)','FontSize',13); xlim([0 40])


legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('C_{HT}, \sigma = 0.1% L_H','FontSize',20);

%% Second iteration: Plot ratio of calculated to known E*
figure
hold on
j = 0; iter = 2;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the very low curvatures 
    subplot(4,3,1)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    subplot(4,3,2)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    subplot(4,3,3)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    subplot(4,3,5)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    subplot(4,3,6)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    subplot(4,3,8)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    subplot(4,3,9)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    
    m = 4; % Now the sharp hilltops 
    subplot(4,3,10)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');
    
    subplot(4,3,11)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');
    
    subplot(4,3,12)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');

end

subplot(4,3,1); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline(1); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline(1); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlim([0 40])
subplot(4,3,5); yline(1); xlim([0 40])
subplot(4,3,6); yline(1); xlim([0 40])
subplot(4,3,7); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlim([0 40])
subplot(4,3,8); yline(1); xlim([0 40])
subplot(4,3,9); yline(1); xlim([0 40])
subplot(4,3,10); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline(1); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline(1); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])

legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('C_{HT} Ratio, \sigma = 0.1% L_H','FontSize',20);

%% Second iteration: Plot standard deviation of E*
figure
hold on
j = 0; iter = 2;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the very low curvatures 
    subplot(4,3,1)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,2)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,3)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,5)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,6)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,8)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,9)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 4; % Now the sharp hilltops 
    subplot(4,3,10)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,11)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,12)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

end

subplot(4,3,1); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline(0); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline(0); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,5); yline(0); xlim([0 40])
subplot(4,3,6); yline(0); xlim([0 40])
subplot(4,3,7); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,8); yline(0); xlim([0 40])
subplot(4,3,9); yline(0); xlim([0 40])
subplot(4,3,10); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline(0); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline(0); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])

legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('Standard Deviation C_{HT}, \sigma = 0.1% L_H','FontSize',20);

%% Third Iteration: A = 5%

iter = 3;  % Which iteration is this?

% Noise Amplitude
Amp = 0.5*Lh;   % Noise Amplitude is 5% of hillslope length 

% White Noise
z_Roer1_R = z_Roer1 + Amp*R;  % Noise matrix
z_Roer10_R = z_Roer10 + Amp*R;
z_Roer30_R = z_Roer30 + Amp*R;
z_Roer100_R = z_Roer100 + Amp*R;


z_Roer1_R = z_Roer1_R + abs(nanmin(z_Roer1_R,[],'all'));    % Reset baselevel as 0
z_Roer10_R = z_Roer10_R + abs(nanmin(z_Roer10_R,[],'all'));    % Reset baselevel as 0
z_Roer30_R = z_Roer30_R + abs(nanmin(z_Roer30_R,[],'all'));       % Reset baselevel as 0
z_Roer100_R = z_Roer100_R + abs(nanmin(z_Roer100_R,[],'all'));       % Reset baselevel as 0


% Pink Noise
z_Roer1_Rp = z_Roer1 + Amp*Rp;
z_Roer10_Rp = z_Roer10 + Amp*Rp;
z_Roer30_Rp = z_Roer30 + Amp*Rp;
z_Roer100_Rp = z_Roer100 + Amp*Rp;


z_Roer1_Rp = z_Roer1_Rp + abs(nanmin(z_Roer1_Rp,[],'all'));    % Reset baselevel as 0
z_Roer10_Rp = z_Roer10_Rp + abs(nanmin(z_Roer10_Rp,[],'all'));    % Reset baselevel as 0
z_Roer30_Rp = z_Roer30_Rp + abs(nanmin(z_Roer30_Rp,[],'all'));       % Reset baselevel as 0
z_Roer100_Rp = z_Roer100_Rp + abs(nanmin(z_Roer100_Rp,[],'all'));       % Reset baselevel as 0


% Red Noise
z_Roer1_Rr = z_Roer1 + Amp*Rr;
z_Roer10_Rr = z_Roer10 + Amp*Rr;
z_Roer30_Rr = z_Roer30 + Amp*Rr;
z_Roer100_Rr = z_Roer100 + Amp*Rr;

z_Roer1_Rr = z_Roer1_Rr + abs(nanmin(z_Roer1_Rr,[],'all'));    % Reset baselevel as 0
z_Roer10_Rr = z_Roer10_Rr + abs(nanmin(z_Roer10_Rr,[],'all'));    % Reset baselevel as 0
z_Roer30_Rr = z_Roer30_Rr + abs(nanmin(z_Roer30_Rr,[],'all'));       % Reset baselevel as 0
z_Roer100_Rr = z_Roer100_Rr + abs(nanmin(z_Roer100_Rr,[],'all'));       % Reset baselevel as 0


j = 0;
for i=5:2:35
    j = j+1;
    a_Las_syn = i/(sqrt(2)*pi*1);
    a_Torr_syn = (i*sqrt(5/2))/(2*pi*1);
    
    % white noise
    [C_Las_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Las_syn,1);
    [C_Torr_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR1] = poly2dgridv4_synth(z_Roer1_R,i);
    
    m = 1; % First row is very low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR1;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR1;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR1;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR1.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR1.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR1.*Lh)./Sc;
    
    [C_Las_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Las_syn,1);
    [C_Torr_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR10] = poly2dgridv4_synth(z_Roer10_R,i);
    
    m = 2; % Second row is low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR10;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR10;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR10;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR10.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR10.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Las_syn,1);
    [C_Torr_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR30] = poly2dgridv4_synth(z_Roer30_R,i);
    
    m = 3; % Third row is sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR30;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR30;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR30;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR30.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR30.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR30.*Lh)./Sc;
   
    
    [C_Las_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Las_syn,1);
    [C_Torr_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR100] = poly2dgridv4_synth(z_Roer100_R,i);
    
    m = 4; % Third row is very sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR100;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR100;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR100;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR100.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR100.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR100.*Lh)./Sc;
    
    % pink noise
    [C_Las_syn_RoerRp1] = conv2_mexh_curv(z_Roer1_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp1] = conv2_mexh_curv(z_Roer1_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp1] = poly2dgridv4_synth(z_Roer1_Rp,i);
    
    m = 1; % First row is very low curvature hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp1;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp1;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp1;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp1.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp1.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp1.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp10] = conv2_mexh_curv(z_Roer10_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp10] = conv2_mexh_curv(z_Roer10_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp10] = poly2dgridv4_synth(z_Roer10_Rp,i);
    
    m = 2; % Second row is gentle hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp10;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp10;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp10;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp10.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp10.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp30] = conv2_mexh_curv(z_Roer30_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp30] = conv2_mexh_curv(z_Roer30_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp30] = poly2dgridv4_synth(z_Roer30_Rp,i);
    
    m = 3; % Third row is sharp hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp30;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp30;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp30;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp30.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp30.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp30.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRp100] = conv2_mexh_curv(z_Roer100_Rp,a_Las_syn,1);
    [C_Torr_syn_RoerRp100] = conv2_mexh_curv(z_Roer100_Rp,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRp100] = poly2dgridv4_synth(z_Roer100_Rp,i);
    
    m = 4; % Third row is sharp hilltops
    C_Las_pink{m,j,iter} = C_Las_syn_RoerRp100;
    C_Torr_pink{m,j,iter} = C_Torr_syn_RoerRp100;
    C_poly_pink{m,j,iter} = C_poly_syn_RoerRp100;
    
    Estar_Las_pink{m,j,iter} = (2*C_Las_syn_RoerRp100.*Lh)./Sc;
    Estar_Torr_pink{m,j,iter} = (2*C_Torr_syn_RoerRp100.*Lh)./Sc;
    Estar_poly_pink{m,j,iter} = (2*C_poly_syn_RoerRp100.*Lh)./Sc;
    
    
    % red noise
    [C_Las_syn_RoerRr1] = conv2_mexh_curv(z_Roer1_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr1] = conv2_mexh_curv(z_Roer1_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr1] = poly2dgridv4_synth(z_Roer1_Rr,i);
    
    m = 1; % First row is moderate curvature hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr1;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr1;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr1;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr1.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr1.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr1.*Lh)./Sc;
    
    [C_Las_syn_RoerRr10] = conv2_mexh_curv(z_Roer10_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr10] = conv2_mexh_curv(z_Roer10_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr10] = poly2dgridv4_synth(z_Roer10_Rr,i);
    
    m = 2; % Second row is gentle hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr10;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr10;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr10;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr10.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr10.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRr30] = conv2_mexh_curv(z_Roer30_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr30] = conv2_mexh_curv(z_Roer30_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr30] = poly2dgridv4_synth(z_Roer30_Rr,i); 
    
    m = 3; % Third row is sharp hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr30;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr30;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr30;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr30.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr30.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr30.*Lh)./Sc;
    
    
    [C_Las_syn_RoerRr100] = conv2_mexh_curv(z_Roer100_Rr,a_Las_syn,1);
    [C_Torr_syn_RoerRr100] = conv2_mexh_curv(z_Roer100_Rr,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerRr100] = poly2dgridv4_synth(z_Roer100_Rr,i); 
    
    m = 4; % Third row is sharp hilltops
    C_Las_red{m,j,iter} = C_Las_syn_RoerRr100;
    C_Torr_red{m,j,iter} = C_Torr_syn_RoerRr100;
    C_poly_red{m,j,iter} = C_poly_syn_RoerRr100;
    
    Estar_Las_red{m,j,iter} = (2*C_Las_syn_RoerRr100.*Lh)./Sc;
    Estar_Torr_red{m,j,iter} = (2*C_Torr_syn_RoerRr100.*Lh)./Sc;
    Estar_poly_red{m,j,iter} = (2*C_poly_syn_RoerRr100.*Lh)./Sc;
    
end

%% Third iteration: Plot Results, High, med, low Cht, change beta
figure
hold on
j = 0; iter = 3;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the moderate curvature 
    subplot(4,3,1)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,2)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,3)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,5)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,6)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,8)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,9)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

    
    m = 4; % Now the very sharp hilltops 
    subplot(4,3,10)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,11)
    hold on
    errorbar(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all'),nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,12)
    hold on
    errorbar(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all'),nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all'),nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

end

subplot(4,3,1); yline((-1)*Cht1,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline((-1)*Cht1,'k--'); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline((-1)*Cht1,'k--'); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline((-1)*Cht10,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,5); yline((-1)*Cht10,'k--'); xlim([0 40])
subplot(4,3,6); yline((-1)*Cht10,'k--'); xlim([0 40])
subplot(4,3,7); yline((-1)*Cht30,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,8); yline((-1)*Cht30,'k--'); xlim([0 40])
subplot(4,3,9); yline((-1)*Cht30,'k--'); xlim([0 40])
subplot(4,3,10); yline((-1)*Cht100,'k--'); ylabel('C_{HT} (m^{-1})','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline((-1)*Cht100,'k--'); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline((-1)*Cht100,'k--'); xlabel('Smoothing scale(m)','FontSize',13); xlim([0 40])



legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('C_{HT}, \sigma = 5% L_H','FontSize',20);

%% Third iteration: Plot ratio of calculated to known E*
figure
hold on
j = 0; iter = 3;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the very low curvatures 
    subplot(4,3,1)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    subplot(4,3,2)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    subplot(4,3,3)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    subplot(4,3,5)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    subplot(4,3,6)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    subplot(4,3,8)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    subplot(4,3,9)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
    
    
    m = 4; % Now the sharp hilltops 
    subplot(4,3,10)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');
    
    subplot(4,3,11)
    hold on
    plot(i,nanmean(C_Las_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_pink{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');
    
    subplot(4,3,12)
    hold on
    plot(i,nanmean(C_Las_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_red{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');

end

subplot(4,3,1); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline(1); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline(1); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlim([0 40])
subplot(4,3,5); yline(1); xlim([0 40])
subplot(4,3,6); yline(1); xlim([0 40])
subplot(4,3,7); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlim([0 40])
subplot(4,3,8); yline(1); xlim([0 40])
subplot(4,3,9); yline(1); xlim([0 40])
subplot(4,3,10); yline(1); ylabel('C_{HT}/known C_{HT}','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline(1); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline(1); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])

legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('C_{HT} Ratio, \sigma = 5% L_H','FontSize',20);

%% Third iteration: Plot standard deviation of E*
figure
hold on
j = 0; iter = 3;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the very low curvatures 
    subplot(4,3,1)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,2)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,3)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,5)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,6)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,8)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,9)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 4; % Now the sharp hilltops 
    subplot(4,3,10)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,11)
    hold on
    plot(i,nanstd(C_Las_pink{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_pink{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_pink{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,12)
    hold on
    plot(i,nanstd(C_Las_red{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_red{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_red{m,j,iter}(:,Lh),0,'all'),'bd');

end

subplot(4,3,1); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); title('White Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,2); yline(0); title('Pink Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,3); yline(0); title('Red Noise Added','FontSize',14); xlim([0 40])
subplot(4,3,4); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,5); yline(0); xlim([0 40])
subplot(4,3,6); yline(0); xlim([0 40])
subplot(4,3,7); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlim([0 40])
subplot(4,3,8); yline(0); xlim([0 40])
subplot(4,3,9); yline(0); xlim([0 40])
subplot(4,3,10); yline(0); ylabel('St. Dev. C_{HT} (m^{-1})','FontSize',13); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,11); yline(0); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])
subplot(4,3,12); yline(0); xlabel('Smoothing scale (m)','FontSize',13); xlim([0 40])

legend('CWT: Lashermes et al.','CWT: Torrence and Compo','Polynomial','FontSize',12); xlim([0 40])
sgtitle('Standard Deviation C_{HT}, \sigma = 5% L_H','FontSize',20);

%% Fourth iteratoin: No noise

iter = 4;  % Which iteration is this?

z_Roer1 = ((K*Sc125^2)/(2*2*eros1))*(log(0.5*(sqrt(1+((2*2*eros1*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros1*xspace_t)/(K*Sc125)).^2)+1);
z_Roer10 = ((K*Sc125^2)/(2*2*eros10))*(log(0.5*(sqrt(1+((2*2*eros10*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros10*xspace_t)/(K*Sc125)).^2)+1);
z_Roer30 = ((K*Sc125^2)/(2*2*eros30))*(log(0.5*(sqrt(1+((2*2*eros30*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros30*xspace_t)/(K*Sc125)).^2)+1);
z_Roer100 = ((K*Sc125^2)/(2*2*eros100))*(log(0.5*(sqrt(1+((2*2*eros100*xspace_t)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros100*xspace_t)/(K*Sc125)).^2)+1);

% Noise (only doing this to keep variable names consistent)
z_Roer1_R = z_Roer1;
z_Roer10_R = z_Roer10;
z_Roer30_R = z_Roer30;
z_Roer100_R = z_Roer100;


z_Roer1_R = z_Roer1_R + abs(nanmin(z_Roer1_R,[],'all'));    % Reset baselevel as 0
z_Roer10_R = z_Roer10_R + abs(nanmin(z_Roer10_R,[],'all'));    % Reset baselevel as 0
z_Roer30_R = z_Roer30_R + abs(nanmin(z_Roer30_R,[],'all'));       % Reset baselevel as 0
z_Roer100_R = z_Roer100_R + abs(nanmin(z_Roer100_R,[],'all'));       % Reset baselevel as 0


j = 0;
for i=5:2:35
    j = j+1;
    a_Las_syn = i/(sqrt(2)*pi*1);
    a_Torr_syn = (i*sqrt(5/2))/(2*pi*1);
    
    %  We will store these data in the white noise category
    [C_Las_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Las_syn,1);
    [C_Torr_syn_RoerR1] = conv2_mexh_curv(z_Roer1_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR1] = poly2dgridv4_synth(z_Roer1_R,i);
    
    m = 1; % First row is very low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR1;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR1;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR1;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR1.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR1.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR1.*Lh)./Sc;
    
    
    [C_Las_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Las_syn,1);
    [C_Torr_syn_RoerR10] = conv2_mexh_curv(z_Roer10_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR10] = poly2dgridv4_synth(z_Roer10_R,i);
    
    m = 2; % Second row is low curvature hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR10;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR10;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR10;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR10.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR10.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR10.*Lh)./Sc;
    
    
    [C_Las_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Las_syn,1);
    [C_Torr_syn_RoerR30] = conv2_mexh_curv(z_Roer30_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR30] = poly2dgridv4_synth(z_Roer30_R,i);
    
    m = 3; % Third row is sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR30;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR30;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR30;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR30.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR30.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR30.*Lh)./Sc;
   
    
    [C_Las_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Las_syn,1);
    [C_Torr_syn_RoerR100] = conv2_mexh_curv(z_Roer100_R,a_Torr_syn,1);
    [~,~,C_poly_syn_RoerR100] = poly2dgridv4_synth(z_Roer100_R,i);
    
    m = 4; % Third row is very sharp hilltops
    C_Las_white{m,j,iter} = C_Las_syn_RoerR100;
    C_Torr_white{m,j,iter} = C_Torr_syn_RoerR100;
    C_poly_white{m,j,iter} = C_poly_syn_RoerR100;
    
    Estar_Las_white{m,j,iter} = (2*C_Las_syn_RoerR100.*Lh)./Sc;
    Estar_Torr_white{m,j,iter} = (2*C_Torr_syn_RoerR100.*Lh)./Sc;
    Estar_poly_white{m,j,iter} = (2*C_poly_syn_RoerR100.*Lh)./Sc;
    
    
end

%% Fourth iteration: Plot Everything
figure
hold on
j = 0; iter = 4;
for i = 5:2:35;
    j = j + 1;
    
    m = 1; % First the moderate curvature 
    subplot(4,3,1)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,2)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht1),'bd');
    
    subplot(4,3,3)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 2; % Gentle hilltops
    subplot(4,3,4)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');

    subplot(4,3,5)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht10),'bd');
    
    subplot(4,3,6)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
   
    
    m = 3; % Now the sharp hilltops 
    subplot(4,3,7)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    subplot(4,3,8)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht30),'bd');
   
    subplot(4,3,9)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
    
    m = 4; % Now the very sharp hilltops 
    subplot(4,3,10)
    hold on
    errorbar(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all'),nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    errorbar(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all'),nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    errorbar(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all'),nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
   
    subplot(4,3,11)
    hold on
    plot(i,nanmean(C_Las_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'ko');
    plot(i,nanmean(C_Torr_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'rs');
    plot(i,nanmean(C_poly_white{m,j,iter}(:,Lh),'all')/((-1)*Cht100),'bd');

    subplot(4,3,12)
    hold on
    plot(i,nanstd(C_Las_white{m,j,iter}(:,Lh),0,'all'),'ko');
    plot(i,nanstd(C_Torr_white{m,j,iter}(:,Lh),0,'all'),'rs');
    plot(i,nanstd(C_poly_white{m,j,iter}(:,Lh),0,'all'),'bd');
    
end

subplot(4,3,1); yline((-1)*Cht1,'k--'); ylabel('C_{HT} (m^{-1}), known E^* = 1'); xlim([0 40]); title('Mean C_{HT} +/- St. Dev');
subplot(4,3,4); yline((-1)*Cht10,'k--'); ylabel('C_{HT} (m^{-1}), known E^* = 10'); xlim([0 40])
subplot(4,3,7); yline((-1)*Cht30,'k--'); ylabel('C_{HT} (m^{-1}), known E^* = 30'); xlim([0 40])
subplot(4,3,10); yline((-1)*Cht100,'k--'); ylabel('C_{HT} (m^{-1}), known E^* = 100'); xlabel('Smoothing scale (m)'); xlim([0 40])


subplot(4,3,2); yline(1); ylabel('C_{HT}/known C_{HT}, E^*=1'); xlim([0 40]); title('Ratio Mean C_{HT}');
subplot(4,3,5); yline(1); ylabel('C_{HT}/known C_{HT}, E^*=10'); xlim([0 40])
subplot(4,3,8); yline(1); ylabel('C_{HT}/known C_{HT}, E^*=30'); xlim([0 40])
subplot(4,3,11); yline(1); ylabel('C_{HT}/known C_{HT}, E^*=100'); xlabel('Smoothing scale (m)'); xlim([0 40])

subplot(4,3,3); yline(0); ylabel('St. Dev. C_{HT} (m^{-1}), E^*=1'); xlim([0 40]); title('Standard Deviation')
subplot(4,3,6); yline(0); ylabel('St. Dev. C_{HT} (m^{-1}), E^*=10'); xlim([0 40])
subplot(4,3,9); yline(0); ylabel('St. Dev. C_{HT} (m^{-1}), E^*=30'); xlim([0 40])
subplot(4,3,12); yline(0); ylabel('St. Dev. C_{HT} (m^{-1}), E^*=100'); xlabel('Smoothing scale (m)'); xlim([0 40])

legend('Lashermes et al.','Torrence and Compo','Polynomial'); xlim([0 40])
sgtitle('No Noise','FontSize',20);

%% Plot Hillslopes for visualization. White Noise
figure
[xplot,yplot] = meshgrid(-100:5:100,-100:5:100,1);  % Make 201x201 matrix (201 so that we have a middle pixel exactly at the hilltop

z_Roer1plot = ((K*Sc125^2)/(2*2*eros1))*(log(0.5*(sqrt(1+((2*2*eros1*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros1*xplot)/(K*Sc125)).^2)+1);
z_Roer10plot = ((K*Sc125^2)/(2*2*eros10))*(log(0.5*(sqrt(1+((2*2*eros10*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros10*xplot)/(K*Sc125)).^2)+1);
z_Roer30plot = ((K*Sc125^2)/(2*2*eros30))*(log(0.5*(sqrt(1+((2*2*eros30*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros30*xplot)/(K*Sc125)).^2)+1);
z_Roer100plot = ((K*Sc125^2)/(2*2*eros100))*(log(0.5*(sqrt(1+((2*2*eros100*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros100*xplot)/(K*Sc125)).^2)+1);

z_Roer1plot = z_Roer1plot + abs(nanmin(z_Roer1plot,[],'all'));    % Reset baselevel as 0
z_Roer10plot = z_Roer10plot+ abs(nanmin(z_Roer10plot,[],'all'));    % Reset baselevel as 0
z_Roer30plot = z_Roer30plot + abs(nanmin(z_Roer30plot,[],'all'));       % Reset baselevel as 0
z_Roer100plot = z_Roer100plot + abs(nanmin(z_Roer100plot,[],'all'));       % Reset baselevel as 0

subplot(5,4,5);
surf(xplot,yplot,z_Roer1plot); title('No Noise Added','FontSize',14)
zlabel('E^*=1');

subplot(5,4,9)
surf(xplot,yplot,z_Roer10plot);
zlabel('E^*=10');

subplot(5,4,13)
surf(xplot,yplot,z_Roer30plot)
zlabel('E^*=30');

subplot(5,4,17)
surf(xplot,yplot,z_Roer100plot)
zlabel('E^*=100');

R_plot = randnd(0,size(xplot,1));

% First 0.1% Lh White Noise
Amp = 0.001*Lh;   % Noise Amplitude is 0.1% of hillslope length 

z_Roer1plot_wpt1 = z_Roer1plot + Amp*R_plot;  % Noise matrix
z_Roer10plot_wpt1 = z_Roer10plot + Amp*R_plot;
z_Roer30plot_wpt1 = z_Roer30plot + Amp*R_plot;
z_Roer100plot_wpt1 = z_Roer100plot + Amp*R_plot;

subplot(5,4,2)
imagesc([-100:5:100],[-100:5:100],Amp*R_plot); col = colorbar; title('White Noise Surface, \sigma = 0.1% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,6);
surf(xplot,yplot,z_Roer1plot_wpt1); title('White Noise, \sigma = 0.1% Lh','FontSize',14)

subplot(5,4,10)
surf(xplot,yplot,z_Roer10plot_wpt1);

subplot(5,4,14)
surf(xplot,yplot,z_Roer30plot_wpt1)

subplot(5,4,18)
surf(xplot,yplot,z_Roer100plot_wpt1)

% 0.5% Lh White Noise
Amp = 0.005*Lh;  % Noise amplitude is 0.5% hillslope length 
z_Roer1plot_w1 = z_Roer1plot + Amp*R_plot;  % Noise matrix
z_Roer10plot_w1 = z_Roer10plot + Amp*R_plot;
z_Roer30plot_w1 = z_Roer30plot + Amp*R_plot;
z_Roer100plot_w1 = z_Roer100plot + Amp*R_plot;

subplot(5,4,3);
imagesc([-100:5:100],[-100:5:100],Amp*R_plot); col = colorbar; title('White Noise Surface, \sigma = 0.5% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,7);
surf(xplot,yplot,z_Roer1plot_w1); title('White Noise, \sigma = 0.5% Lh','FontSize',14)

subplot(5,4,11)
surf(xplot,yplot,z_Roer10plot_w1);

subplot(5,4,15)
surf(xplot,yplot,z_Roer30plot_w1)

subplot(5,4,19)
surf(xplot,yplot,z_Roer100plot_w1)


% 5% Lh White Noise
Amp = 0.05*Lh;  % Noise amplitude is 0.5% hillslope length 
z_Roer1plot_w10 = z_Roer1plot + Amp*R_plot;  % Noise matrix
z_Roer10plot_w10 = z_Roer10plot + Amp*R_plot;
z_Roer30plot_w10 = z_Roer30plot + Amp*R_plot;
z_Roer100plot_w10 = z_Roer100plot + Amp*R_plot;

subplot(5,4,4);
imagesc([-100:5:100],[-100:5:100],Amp*R_plot); col = colorbar; title('White Noise Surface, \sigma = 5% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,8);
surf(xplot,yplot,z_Roer1plot_w10); title('White Noise, \sigma = 5% Lh','FontSize',14)

subplot(5,4,12)
surf(xplot,yplot,z_Roer10plot_w10);

subplot(5,4,16)
surf(xplot,yplot,z_Roer30plot_w10)

subplot(5,4,20)
surf(xplot,yplot,z_Roer100plot_w10)

%% Plot Hillslopes for visualization. Pink Noise
figure
[xplot,yplot] = meshgrid(-100:5:100,-100:5:100,1);  % Make 201x201 matrix (201 so that we have a middle pixel exactly at the hilltop

z_Roer1plot = ((K*Sc125^2)/(2*2*eros1))*(log(0.5*(sqrt(1+((2*2*eros1*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros1*xplot)/(K*Sc125)).^2)+1);
z_Roer10plot = ((K*Sc125^2)/(2*2*eros10))*(log(0.5*(sqrt(1+((2*2*eros10*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros10*xplot)/(K*Sc125)).^2)+1);
z_Roer30plot = ((K*Sc125^2)/(2*2*eros30))*(log(0.5*(sqrt(1+((2*2*eros30*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros30*xplot)/(K*Sc125)).^2)+1);
z_Roer100plot = ((K*Sc125^2)/(2*2*eros100))*(log(0.5*(sqrt(1+((2*2*eros100*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros100*xplot)/(K*Sc125)).^2)+1);

z_Roer1plot = z_Roer1plot + abs(nanmin(z_Roer1plot,[],'all'));    % Reset baselevel as 0
z_Roer10plot = z_Roer10plot+ abs(nanmin(z_Roer10plot,[],'all'));    % Reset baselevel as 0
z_Roer30plot = z_Roer30plot + abs(nanmin(z_Roer30plot,[],'all'));       % Reset baselevel as 0
z_Roer100plot = z_Roer100plot + abs(nanmin(z_Roer100plot,[],'all'));       % Reset baselevel as 0

subplot(5,4,5);
surf(xplot,yplot,z_Roer1plot); title('No Noise Added','FontSize',14)
zlabel('E^*=1');

subplot(5,4,9)
surf(xplot,yplot,z_Roer10plot);
zlabel('E^*=10');

subplot(5,4,13)
surf(xplot,yplot,z_Roer30plot)
zlabel('E^*=30');

subplot(5,4,17)
surf(xplot,yplot,z_Roer100plot)
zlabel('E^*=100');

Rp_plot = randnd(-1,size(xplot,1));

% First 0.1% Lh Pink Noise
Amp = 0.001*Lh;   % Noise Amplitude is 0.1% of hillslope length 

z_Roer1plot_ppt1 = z_Roer1plot + Amp*Rp_plot;  % Noise matrix
z_Roer10plot_ppt1 = z_Roer10plot + Amp*Rp_plot;
z_Roer30plot_ppt1 = z_Roer30plot + Amp*Rp_plot;
z_Roer100plot_ppt1 = z_Roer100plot + Amp*Rp_plot;

subplot(5,4,2)
imagesc([-100:5:100],[-100:5:100],Amp*Rp_plot); col = colorbar; title('Pink Noise Surface, \sigma = 0.1% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,6);
surf(xplot,yplot,z_Roer1plot_ppt1); title('Pink Noise, \sigma = 0.1% Lh','FontSize',14)

subplot(5,4,10)
surf(xplot,yplot,z_Roer10plot_ppt1);

subplot(5,4,14)
surf(xplot,yplot,z_Roer30plot_ppt1)

subplot(5,4,18)
surf(xplot,yplot,z_Roer100plot_ppt1)

% 0.5% Lh White Noise
Amp = 0.005*Lh;  % Noise amplitude is 0.5% hillslope length 
z_Roer1plot_p1 = z_Roer1plot + Amp*Rp_plot;  % Noise matrix
z_Roer10plot_p1 = z_Roer10plot + Amp*Rp_plot;
z_Roer30plot_p1 = z_Roer30plot + Amp*Rp_plot;
z_Roer100plot_p1 = z_Roer100plot + Amp*Rp_plot;

subplot(5,4,3);
imagesc([-100:5:100],[-100:5:100],Amp*Rp_plot); col = colorbar; title('Pink Noise Surface, \sigma = 0.5% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,7);
surf(xplot,yplot,z_Roer1plot_p1); title('Pink Noise, \sigma = 0.5% Lh','FontSize',14)

subplot(5,4,11)
surf(xplot,yplot,z_Roer10plot_p1);

subplot(5,4,15)
surf(xplot,yplot,z_Roer30plot_p1)

subplot(5,4,19)
surf(xplot,yplot,z_Roer100plot_p1)


% 5% Lh Pink Noise
Amp = 0.05*Lh;  % Noise amplitude is 5% hillslope length 
z_Roer1plot_p10 = z_Roer1plot + Amp*Rp_plot;  % Noise matrix
z_Roer10plot_p10 = z_Roer10plot + Amp*Rp_plot;
z_Roer30plot_p10 = z_Roer30plot + Amp*Rp_plot;
z_Roer100plot_p10 = z_Roer100plot + Amp*Rp_plot;

subplot(5,4,4);
imagesc([-100:5:100],[-100:5:100],Amp*Rp_plot); col = colorbar; title('Pink Noise Surface, \sigma = 5% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,8);
surf(xplot,yplot,z_Roer1plot_p10); title('Pink Noise, \sigma = 5% Lh','FontSize',14)

subplot(5,4,12)
surf(xplot,yplot,z_Roer10plot_p10);

subplot(5,4,16)
surf(xplot,yplot,z_Roer30plot_p10)

subplot(5,4,20)
surf(xplot,yplot,z_Roer100plot_p10)

%% Plot Hillslopes for visualization. Red Noise
figure
[xplot,yplot] = meshgrid(-100:5:100,-100:5:100,1);  % Make 201x201 matrix (201 so that we have a middle pixel exactly at the hilltop

z_Roer1plot = ((K*Sc125^2)/(2*2*eros1))*(log(0.5*(sqrt(1+((2*2*eros1*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros1*xplot)/(K*Sc125)).^2)+1);
z_Roer10plot = ((K*Sc125^2)/(2*2*eros10))*(log(0.5*(sqrt(1+((2*2*eros10*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros10*xplot)/(K*Sc125)).^2)+1);
z_Roer30plot = ((K*Sc125^2)/(2*2*eros30))*(log(0.5*(sqrt(1+((2*2*eros30*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros30*xplot)/(K*Sc125)).^2)+1);
z_Roer100plot = ((K*Sc125^2)/(2*2*eros100))*(log(0.5*(sqrt(1+((2*2*eros100*xplot)/(K*Sc125)).^2)+1))...
    -sqrt(1+((2*2*eros100*xplot)/(K*Sc125)).^2)+1);

z_Roer1plot = z_Roer1plot + abs(nanmin(z_Roer1plot,[],'all'));    % Reset baselevel as 0
z_Roer10plot = z_Roer10plot+ abs(nanmin(z_Roer10plot,[],'all'));    % Reset baselevel as 0
z_Roer30plot = z_Roer30plot + abs(nanmin(z_Roer30plot,[],'all'));       % Reset baselevel as 0
z_Roer100plot = z_Roer100plot + abs(nanmin(z_Roer100plot,[],'all'));       % Reset baselevel as 0

subplot(5,4,5);
surf(xplot,yplot,z_Roer1plot); title('No Noise Added','FontSize',14)
zlabel('E^*=1');

subplot(5,4,9)
surf(xplot,yplot,z_Roer10plot);
zlabel('E^*=10');

subplot(5,4,13)
surf(xplot,yplot,z_Roer30plot)
zlabel('E^*=30');

subplot(5,4,17)
surf(xplot,yplot,z_Roer100plot)
zlabel('E^*=100');

Rr_plot = randnd(-2,size(xplot,1));

% First 0.1% Lh Red Noise
Amp = 0.001*Lh;   % Noise Amplitude is 0.1% of hillslope length 

z_Roer1plot_rpt1 = z_Roer1plot + Amp*Rr_plot;  % Noise matrix
z_Roer10plot_rpt1 = z_Roer10plot + Amp*Rr_plot;
z_Roer30plot_rpt1 = z_Roer30plot + Amp*Rr_plot;
z_Roer100plot_rpt1 = z_Roer100plot + Amp*Rr_plot;

subplot(5,4,2)
imagesc([-100:5:100],[-100:5:100],Amp*Rr_plot); col = colorbar; title('Red Noise Surface, \sigma = 0.1% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,6);
surf(xplot,yplot,z_Roer1plot_rpt1); title('Red Noise, \sigma = 0.1% Lh','FontSize',14)

subplot(5,4,10)
surf(xplot,yplot,z_Roer10plot_rpt1);

subplot(5,4,14)
surf(xplot,yplot,z_Roer30plot_rpt1)

subplot(5,4,18)
surf(xplot,yplot,z_Roer100plot_rpt1)

%0.5% Lh Red Noise
Amp = 0.005*Lh;  % Noise amplitude is 0.5% hillslope length 
z_Roer1plot_r1 = z_Roer1plot + Amp*Rr_plot;  % Noise matrix
z_Roer10plot_r1 = z_Roer10plot + Amp*Rr_plot;
z_Roer30plot_r1 = z_Roer30plot + Amp*Rr_plot;
z_Roer100plot_r1 = z_Roer100plot + Amp*Rr_plot;

subplot(5,4,3);
imagesc([-100:5:100],[-100:5:100],Amp*Rr_plot); col = colorbar; title('Red Noise Surface, \sigma = 0.5% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,7);
surf(xplot,yplot,z_Roer1plot_r1); title('Red Noise, \sigma = 0.5% Lh','FontSize',14)

subplot(5,4,11)
surf(xplot,yplot,z_Roer10plot_r1);

subplot(5,4,15)
surf(xplot,yplot,z_Roer30plot_r1)

subplot(5,4,19)
surf(xplot,yplot,z_Roer100plot_r1)


% 5% Lh Red Noise
Amp = 0.05*Lh;  % Noise amplitude is 5% hillslope length 
z_Roer1plot_r10 = z_Roer1plot + Amp*Rr_plot;  % Noise matrix
z_Roer10plot_r10 = z_Roer10plot + Amp*Rr_plot;
z_Roer30plot_r10 = z_Roer30plot + Amp*Rr_plot;
z_Roer100plot_r10 = z_Roer100plot + Amp*Rr_plot;

subplot(5,4,4);
imagesc([-100:5:100],[-100:5:100],Amp*Rr_plot); col = colorbar; title('Red Noise Surface, \sigma = 5% Lh','FontSize',14); axis image;
ylabel(col,'Noise Amplitude'); set(col,'FontSize',12);

subplot(5,4,8);
surf(xplot,yplot,z_Roer1plot_r10); title('Red Noise, \sigma = 5% Lh','FontSize',14)

subplot(5,4,12)
surf(xplot,yplot,z_Roer10plot_r10);

subplot(5,4,16)
surf(xplot,yplot,z_Roer30plot_r10)

subplot(5,4,20)
surf(xplot,yplot,z_Roer100plot_r10)

