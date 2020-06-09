function [xx,yy,zz]=ellipsoid3(varargin)
%ELLIPSOID3 Generate 3D rotated ellipsoid.
%   [X,Y,Z]=ELLIPSOID3(AX,AY,AZ,RX,RY,RZ,N) generates three
%   (N+1)-by-(N+1) matrices so that SURF(X,Y,Z) produces an
%   ellipsoid centered at (0,0,0) with semi-axis radii AX, AY, AZ, 
%	and rotation angles RX, RY, RZ around each axis (in degree).
% 
%   [X,Y,Z]=ELLIPSOID3(AX,AY,AZ,RX,RY,RZ) uses N = 20.
%
%   ELLIPSOID3(...) and ELLIPSOID3(...,N) with no output arguments
%   graph the ellipsoid as a SURFACE and do not return anything.
%
%   ELLIPSOID3(AX,...) plots into AX instead of GCA.
%
%   The ellipsoidal data is generated using the equation:
%
%       X^2          Y^2          Z^2
%    --------  +  --------  +  --------  =  1
%      AX^2         AY^2         AY^2
%
%	thus applies a 3D rotation using RX, RY and RZ angles.
%
%   See also SPHERE, ELLIPSOID, CYLINDER.
%
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%	Created: 2020-06-07 in Yogyakarta, Indonesia
%
%	Based on the original script by Laurens Schalekamp and Damian T. Packer
%   Copyright 1984-2002 The MathWorks, Inc. 

% Parse possible Axes input
narginchk(6,8);
[cax,args,nargs] = axescheck(varargin{:});

[ax,ay,az,rx,ry,rz] = deal(args{1:6});
n  = 20;

if nargs > 6
	n = args{7}; 
end

% generates a unitary sphere
[x,y,z] = sphere(n);

% extends semi-axis radii
x = ax*x;
y = ay*y;
z = az*z;

% rotation matrix
Rx = [1 0 0;0 cosd(rx) sind(rx);0 -sind(rx) cosd(rx)];
Ry = [cosd(ry) 0 -sind(ry);0 1 0;sind(ry) 0 cosd(ry)];
Rz = [cosd(rz) sind(rz) 0;-sind(rz) cosd(rz) 0;0 0 1];
R = Rz*Ry*Rx;

xyz = R*[x(:),y(:),z(:)]';

sz = size(x);
x = reshape(xyz(1,:),sz);
y = reshape(xyz(2,:),sz);
z = reshape(xyz(3,:),sz);

if(nargout == 0)
    cax = newplot(cax);
	surf(x,y,z,'parent',cax)
else
	xx = x;
	yy = y;
	zz = z;
end
