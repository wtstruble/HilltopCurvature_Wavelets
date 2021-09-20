function [ELEV,GRAD,CURV]=poly2dgridv4(DEM,w)

%This calculates gradient and curvature for interior nodes in a dem using
%the polynomial:  z = ax^2 + by^2 + cxy + dx + ey + f.  The subroutine
%makenormeqn.m solves the local matrix inversion with a weighting least
%squares algorithm (Jo Wood's PhD thesis, 1996, online). 
%w is the window size (in nodes...must be odd number!)

%13-Aug-2010  
%2nd version takes the Arc .flt  files as input and then pumps 
%out arc .flt raster ouputs (elev, grad, curv) with the same header and
%no data values.  

%5-jul-2011
%changed input/output to accept .flt files as per adam booth's ".flt" 
%scripts

%9-Aug-2018
% Will Struble changed the format of the code to be compatible with
% TopoToolbox. Code now maintains spatial reference of the original DEM,
% and improves upon the gradient and curvature functions of TopoToolbox
% (Schwangart). 

% Spring 2021
% Will Struble modified code to be compatible with synethic hillslope code
% (in practice, similar to v. 3).

tic     %start timer

%set cellsize
%g=DEM.cellsize;
g = 1;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this creates the left side of the normal equations for matrix solutions
[nmat,X2,Y2,X,Y,weight]=makenormeqn(w,g);

%This calculates how many nodes from the edge of the matrix you can start 
%doing this fit...
edge=(w-1)/2;

%calculate the dem size
%[m n]=size(DEM.Z);
[m n] = size(DEM);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%assign NaN values to the fringe of the null grid and 1's to interior
fringevalue=NaN;
null(1:m,1:n)=fringevalue;
for i=(edge+1):(m-edge)
    for j=(edge+1):(n-edge)
        null(i,j)=1;
    end
end

%create the output grids and insert "0" values...
elev(1:m,1:n)=0;grad(1:m,1:n)=0;curv(1:m,1:n)=0;

%loop through the big grid and do the local fit
for i=(edge+1):(m-edge)
    for j=(edge+1):(n-edge)
        %z=DEM.Z(i-edge:i+edge,j-edge:j+edge);    %%%%%%%%%%%%%%%%%%%%%%       
        z = DEM(i-edge:i+edge,j-edge:j+edge);
        %Z components of the normal equations
        zvec(1:6)=0;
        zvec(1)=sum(sum(z.*weight.*X2));
        zvec(2)=sum(sum(z.*weight.*Y2));
        zvec(3)=sum(sum(z.*weight.*X.*Y));
        zvec(4)=sum(sum(z.*weight.*X));
        zvec(5)=sum(sum(z.*weight.*Y));
        zvec(6)=sum(sum(z.*weight));
                
        %Do the matrix inversion
        A=nmat\zvec';
        
        %calculate elev, gradient, and curvature based on derivatives
        elev(i,j)=A(6);
        grad(i,j)=sqrt(A(4).^2+A(5).^2);
        curv(i,j)=2*(A(1)+A(2));
    end
end

%add null values to the fringes of these grids and then change NaN values
%to NODATA_value at all points in the matrix...
elev=elev.*null;
elev(isnan(elev))=NaN;
grad=grad.*null;
grad(isnan(grad))=NaN;
curv=curv.*null;
curv(isnan(curv))=NaN;

% ELEV=DEM;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRAD=DEM;
% CURV=DEM;
CURV = curv;
% ELEV.Z=elev;
% GRAD.Z=grad;
% CURV.Z=curv;
ELEV = elev;

toc     %stop timer