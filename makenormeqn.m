function [nmat,X2,Y2,X,Y,weight]=makenormeqn(win,g)

%This calculates the left side normal equation for fitting a 2-d 2nd order
%polynomial to gridded DEM data.  It includes a weighting function.

% z is the local dem matrix with nxn dimensions
% A is the vector with coefficients for the following polynomial fit:
% z = ax^2 + by^2 + cxy + dx + ey + f;

%This calculates how far (distance units) from the edge of the matrix you 
%can start doing this fit...nodes*grid spacing
edge=g*(win-1)/2;

%This creates a mesh of x,y values in the coordinates based on the window
%size chosen
x=(-1*edge):g:edge;
y=(-1*edge):g:edge;
[X,Y]=meshgrid(x,y);
Y=flipud(Y);
X2=X.^2; Y2=Y.^2; X3=X.^3; Y3=Y.^3; X4=X.^4; Y4=Y.^4;

%Define the weighting matrix according to a nonlinear function
distance=sqrt((X.^2)+(Y.^2));
dcrit=sqrt(edge^2+edge^2);
weight=1./((distance./dcrit)+1).^(3/2);

%make the normal equations and matrix
%fill left side of normal equation matrix
nmat(1:6,1:6)=0;

%set-up the unique matrix components and then fill it in...
x4w=sum(sum(X4.*weight));
y4w=sum(sum(Y4.*weight));
x2y2w=sum(sum(X2.*Y2.*weight));
x2w=sum(sum(X2.*weight));
y2w=sum(sum(Y2.*weight));
ww=sum(sum(weight));
x3yw=sum(sum(X3.*Y.*weight));
x3w=sum(sum(X3.*weight));
x2yw=sum(sum(X2.*Y.*weight));
xy3w=sum(sum(X.*Y3.*weight));
xy2w=sum(sum(X.*Y2.*weight));
y3w=sum(sum(Y3.*weight));
xyw=sum(sum(X.*Y.*weight));
xw=sum(sum(X.*weight));
yw=sum(sum(Y.*weight));

%fill in the matrix...
nmat(1,1)=x4w;
nmat(1,2)=x2y2w;
nmat(1,3)=x3yw;
nmat(1,4)=x3w;
nmat(1,5)=x2yw;
nmat(1,6)=x2w;

nmat(2,1)=x2y2w;
nmat(2,2)=y4w;
nmat(2,3)=xy3w;
nmat(2,4)=xy2w;
nmat(2,5)=y3w;
nmat(2,6)=y2w;

nmat(3,1)=x3yw;
nmat(3,2)=xy3w;
nmat(3,3)=x2y2w;
nmat(3,4)=x2yw;
nmat(3,5)=xy2w;
nmat(3,6)=xyw;

nmat(4,1)=x3w;
nmat(4,2)=xy2w;
nmat(4,3)=x2yw;
nmat(4,4)=x2w;
nmat(4,5)=xyw;
nmat(4,6)=xw;

nmat(5,1)=x2yw;
nmat(5,2)=y3w;
nmat(5,3)=xy2w;
nmat(5,4)=xyw;
nmat(5,5)=y2w;
nmat(5,6)=yw;

nmat(6,1)=x2w;
nmat(6,2)=y2w;
nmat(6,3)=xyw;
nmat(6,4)=xw;
nmat(6,5)=yw;
nmat(6,6)=ww;










