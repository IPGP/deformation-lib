function varargout=pcdmv(varargin)
%function [ue,un,uv]=pcdmv(X,Y,depth,omegaX,omegaY,omegaZ,DVtot,A,B,nu)
%PCDMV	Point source compound dislocation model.
%	[uE,uN,uV]=PCDMV(X,Y,DEPTH,OMEGAX,OMEGAY,OMEGAZ,DVTOT,A,B,NU) calculates
%	the surface displacements associated with a point Compound Dislocation 
%	Model by Nikkhoo (2017) that is composed of three mutually orthogonal 
%	point tensile dislocations in an elastic half-space.
%
%	This model is able to approximate displacement due to any shape of 
%	source in the far field: isotropic, sill, dyke, pipe, or any ellipsoid
%	in all orientations in space.
%
%	Equations have been fully vectorized so all input parameter can be
%	scalars, vectors or N-D matrix of the same number of elements. Output 
%	parameters will have the same size as X.
%	
%
%	--- Input parameters
%
%	X and Y:
%	   Horizontal coordinates (East, North) of calculation points relative
%	   to source located at (0,0).
%
%	DEPTH:
%	   Depth of the source from calculation points, same unit as X and Y.
%	   Note that you might add source depth to elevation at each 
%	   calculation points to approximate the topographic effects.
%
%	OMEGAX, OMEGAY and OMEGAZ:
%	   Clockwise rotation angles about X, Y and Z axes, respectively, that 
%	   specify the orientation of the pCDM in space, in degrees.
%
%	DVTOT:
%	   Total potency DVTOT = DVX+DVY+DVZ of the PTDs that before applying 
%	   the rotations are normal to the X, Y and Z axes, respectively. The 
%	   potency has the unit of volume (the unit of displacements and CDM 
%	   semi-axes to the power of 3).
%
%	A:
%	   Horizontal over total volume variation ratio A = DVZ/DVTOT
%
%	B:
%	   Vertical volume variation ratio B = DVY/(DVX+DVY)
%
%	   Examples:
% 		A = 1,      any B : horizontal sill
% 		A = 0,    B = 0.5 : vertical pipe
% 		A = 0, B = 0 or 1 : vertical dyke
% 		A = 1/3,  B = 0.5 : isotropic source
% 		A > 0,      B > 0 : any ellipsoid
%
%	NU:
%	   Poisson's ratio, optional and dimensionless (default is 0.25 for 
%	   an isotropic medium).
%	
%	
%	--- Outputs parameters
%
%	uE, uN and uV:
%	   Calculated displacement vector components in EFCS. Will have the 
%	   same unit as X, Y and DEPTH in inputs.
%
%
%	--- Glossary
%
%	   pCDM: point Compound Dislocation Model
%	   PTD: Point Tensile Dislocation
%	   EFCS: Earth-Fixed Coordinate System
%
%
%	--- References
%
%	   Nikkhoo, M., Walter, T. R., Lundgren, P. R., Prats-Iraola, P. (2017):
%	    Compound dislocation models (CDMs) for volcano deformation analyses.
%	    Geophys. J. Int., 208(2): 877-894, doi:10.1093/gji/ggw427
%
%	   http://www.volcanodeformation.com
%
%
%	Authors: François Beauducel, Antoine Villié, and Mehdi Nikkhoo
%	Created: 2015-05-22 in GFZ Potsdam (Germany) by Mehdi Nikkhoo
%	Modified: 2018-07-18 in UGM Yogyakarta (Indonesia) by Antoine Villié
%	Updated: 2019-08-21

%	Copyright (c) 2016 Mehdi Nikkhoo
%	Copyright (c) 2018 Antoine Villié
%	Copyright (c) 2019 François Beauducel
%
%	Permission is hereby granted, free of charge, to any person obtaining a
%	copy of this software and associated documentation files
%	(the "Software"), to deal in the Software without restriction, including
%	without limitation the rights to use, copy, modify, merge, publish,
%	distribute, sublicense, and/or sell copies of the Software, and to permit
%	persons to whom the Software is furnished to do so, subject to the
%	following conditions:
%
%	The above copyright notice and this permission notice shall be included
%	in all copies or substantial portions of the Software.
%
%	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%	OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
%	NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
%	DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
%	OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
%	USE OR OTHER DEALINGS IN THE SOFTWARE.

% checks input arguments
if length(find(unique(cellfun(@numel,varargin)) ~= 1)) > 1
	error('All input arguments must be scalars or have the same numel or elements.')
end
sx = size(varargin{1});

% inputs arguments are converted to column vectors
X  = varargin{1}(:);
Y  = varargin{2}(:);
D  = varargin{3}(:);
OX = varargin{4}(:);
OY = varargin{5}(:);
OZ = varargin{6}(:);
DV = varargin{7}(:);
A  = varargin{8}(:);
B  = varargin{9}(:);
if nargin < 10
	NU = 0.25;
else
	NU = varargin{10}(:);
end

% recomputes DV per component (original Nikkhoo's convention) from DVTOT, A and B parameters
DVz = DV.*A;
DVy = (DV-DVz).*B;
DVx = (DV-DVz).*(1-B);

% R1, R2 and R3 are the coefficients from the 3-D matrix of rotation

% calculates contribution of the first PTD
R1 =  cosd(OY).*cosd(OZ);
R2 = -cosd(OY).*sind(OZ);
R3 =  sind(OY);
nstrike = sqrt(R1.^2 + R2.^2);
strike = atan2d(-R2./nstrike,R1./nstrike);
strike(isnan(strike)) = 0;
[e1,n1,v1] = PTD(X,Y,D,strike-90,acosd(R3),DVx,NU);

% calculates contribution of the second PTD
R1 =  cosd(OZ).*sind(OX).*sind(OY) + cosd(OX).*sind(OZ);
R2 = -sind(OX).*sind(OY).*sind(OZ) + cosd(OX).*cosd(OZ);
R3 = -sind(OX).*cosd(OY);
nstrike = sqrt(R1.^2 + R2.^2);
strike = atan2d(-R2./nstrike,R1./nstrike);
strike(isnan(strike)) = 0;
[e2,n2,v2] = PTD(X,Y,D,strike-90,acosd(R3),DVy,NU);

% calculates contribution of the third PTD
R1 = -cosd(OX).*sind(OY).*cosd(OZ) + sind(OX).*sind(OZ);
R2 =  cosd(OX).*sind(OY).*sind(OZ) + sind(OX).*cosd(OZ);
R3 =  cosd(OX).*cosd(OY);
nstrike = sqrt(R1.^2 + R2.^2);
strike = atan2d(-R2./nstrike,R1./nstrike);
strike(isnan(strike)) = 0;
[e3,n3,v3] = PTD(X,Y,D,strike-90,acosd(R3),DVz,NU);

% reshapes outputs to the original size of X
varargout{1} = reshape(e1+e2+e3,sx);
varargout{2} = reshape(n1+n2+n3,sx);
varargout{3} = reshape(v1+v2+v3,sx);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ue,un,uv]=PTD(x,y,d,beta,dip,dv,nu)
% PTD calculates surface displacements associated with a tensile 
% point dislocation (PTD) in an elastic half-space (Okada, 1985).

tmp = x.*cosd(beta) - y.*sind(beta);
y =   x.*sind(beta) + y.*cosd(beta);
x = tmp;
r = sqrt(x.^2 + y.^2 + d.^2);
q25 = 3*(y.*sind(dip) - d.*cosd(dip)).^2./r.^5;

e = x.*(q25 - sind(dip).^2.*(1-2*nu).*(1./r.^3 - 1./(r.*(r+d).^2) + y.^2.*(3*r+d)./(r.^3.*(r+d).^3)));
n = y.*(q25 - sind(dip).^2.*(1-2*nu).*(1./(r.*(r+d).^2)           - x.^2.*(3*r+d)./(r.^3.*(r+d).^3)));
v = (d.*q25 - sind(dip).^2.*(1-2*nu).*(1./(r.*(r+d))              - x.^2.*(2*r+d)./(r.^3.*(r+d).^2)));

ue = ( e.*cosd(beta)+n.*sind(beta)).*dv/2/pi;
un = (-e.*sind(beta)+n.*cosd(beta)).*dv/2/pi;
uv = v.*dv/2/pi;

end
