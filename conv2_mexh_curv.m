function [C,frq,wave] = conv2_mexh_curv(dem,a,dx)

%[C,frq,wave] = conv2_mexh(dem,a,dx);
%Computes the 2D CWT of dem using the Mexican Hat wavelet.
%
%dem = digital elevation model
%a = wavelet scale
%dx = grid spacing (same in x- and y-directions)
%
%C = array of wavelet coefficients, indexed to dem
%frq = bandpass frequency of wavelet at scale a
%wave = wavelength (inverse of frq)
%
% Written by A.M. Booth, Portland State.
% Derived from Automated Landslide Mapping Toolkit
% (http://web.pdx.edu/~boothad/tools.html)

tic     %start timer
            
%Generate the mexican hat wavelet kernel at wavelet scale a.  The kernel 
%must be large enough for the wavelet to decay to ~0 at the edges.
[X,Y] = meshgrid(-8*a:8*a,-8*a:8*a);
%This psi has not been scaled, so it is just proportional to curvature:
%psi = (1/a).*(2 - (X/a).^2 - (Y/a).^2).*exp(-((X/a).^2 + (Y/a).^2)/2);
%This psi has been scaled, so it is equal to curvature:
psi = (-1/(pi*(a*dx)^4))*(1 - (X.*X + Y.*Y)/(2*a^2)).*exp(-(X.*X + Y.*Y)/(2*a^2));

%Convolve dem with psi using matlab's conv2 function, multiplying by dx^2
%to approximate the double integral.  'same' crops C to same size as dem.
C = (dx^2)*conv2(dem,psi,'same');

%Mask edge effects with NaN (no data) values.
[nrows,ncols] = size(C);
fringeval = ceil(a*4);
C(1:fringeval,:) = NaN;
C(:,1:fringeval) = NaN;
C(nrows-fringeval+1:nrows,:) = NaN;
C(:,ncols-fringeval+1:ncols) = NaN;

%Frequency and wavelength vectors:
wave = 2*pi*dx*a/(5/2)^(1/2);      %Torrence and Compo [1998]
frq = 1./wave;

toc     %stop timer