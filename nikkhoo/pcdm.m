function [ue,un,uv]=pcdm(X,Y,depth,omegaX,omegaY,omegaZ,DVtot,A,B,nu)
%PCDM	Point source compound dislocation model.
%	PCDM(X,Y,DEPTH,OMEGAX,OMEGAY,OMEGAZ,DVTOT,A,B,NU) calculates the 
%	surface displacements associated with a point Compound Dislocation 
%	Model by Nikkhoo (2016) that is composed of three mutually orthogonal 
%	point tensile dislocations in an elastic half-space.
%
%	This model is able to approximate displacement due to any shape of 
%	source in the far field: isotropic, sill, dyke, pipe, or any ellipsoid
%	in all orientations in space.
%	
%	--- Input parameters
%
%	X and Y:
%	   Horizontal coordinates (East, North) of calculation points relative
%	   to source located at (0,0). X and Y must have the same size.
%
%	DEPTH:
%	   Depth of the source from calculation points, scalar of same size and
%	   same unit as X and Y. Note that you might add source depth to elevation
%	   at each calculation points to approximate the topographic effects.
%
%	omegaX, omegaY and omegaZ:
%	   Clockwise rotation angles about X, Y and Z axes, respectively, that 
%	   specify the orientation of the pCDM in space. The input values must
%	   be scalar expressed in degrees.
%
%	DVtot:
%	   Total potency DVtot = DVx+DVy+DVz of the PTDs that before applying 
%	   the rotations are normal to the X, Y and Z axes, respectively. The 
%	   potency has the unit of volume (the unit of displacements and CDM 
%	   semi-axes to the power of 3).
%
%	A:
%	   Horizontal over total volume variation ratio A = dVZ/(dVX+dVY+dVZ)
%
%	B:
%	   Vertical volume variation ratio B = dVY/(dVX+dVY)
%
%	   Examples:
% 		A = 1, any B : horizontal sill
% 		A = 0, B = 0.5 : vertical pipe
% 		A = 0, B = 0 or 1 : vertical dyke
% 		A = 1/3, B = 0.5 : isotropic source
% 		A > 0, B > 0 : any ellipsoid
%
%	nu:
%	   Poisson's ratio (default is 0.25 for isotropic medium).
%	
%	
%	--- Outputs parameters
%
%	ue, un and uv:
%	   Calculated displacement vector components in EFCS. ue, un and uv have
%	   the same unit as X, Y and DEPTH in inputs.
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
%	   Nikkhoo, M., Walter, T. R., Lundgren, P. R., Prats-Iraola, P. (2016):
%	    Compound dislocation models (CDMs) for volcano deformation analyses.
%	    Geophys. J. Int.
%
%	   http://www.volcanodeformation.com
%
%
%	Authors: Mehdi Nikkhoo, modified by Antoine Villié and François Beauducel
%	Created: 2015-05-22 in GFZ Potsdam (Germany)
%	Modified: 2018-07-18 in UGM Yogyakarta (Indonesia)
%	Updated: 2019-08-05

%	Copyright (c) 2016 Mehdi Nikkhoo
%	Copyright (c) 2018 Antoine Villié and François Beauducel
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
if nargin < 10
	nu = 0.25;
end
if isscalar(depth)
	depth = repmat(depth,size(X));
end
if length(unique([numel(X),numel(Y),numel(depth)])) > 1
	error('X, Y and depth must have the same size.')
end
if ~all([isscalar(omegaX),isscalar(omegaY),isscalar(omegaZ),isscalar(A),isscalar(B),isscalar(nu)])
	error('OMEGAX, OMEGAY, OMEGAZ, DVTOT, A, B, and NU must be scalars.')
end

% recomputes DV per component (original Nikkhoo's convention) from DVtot, A and B parameters
DVz = DVtot*A;
DVy = (DVtot-DVz)*B;
DVx = (DVtot-DVz)*(1-B);

X = X(:);
Y = Y(:);
depth = depth(:);

Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

Vstrike1 = [-R(2,1),R(1,1),0];
Vstrike1 = Vstrike1/norm(Vstrike1);
strike1 = atan2(Vstrike1(1),Vstrike1(2))*180/pi;
if isnan(strike1)
    strike1 = 0;
end
dip1 = acosd(R(3,1));

Vstrike2 = [-R(2,2),R(1,2),0];
Vstrike2 = Vstrike2/norm(Vstrike2);
strike2 = atan2(Vstrike2(1),Vstrike2(2))*180/pi;
if isnan(strike2)
    strike2 = 0;
end
dip2 = acosd(R(3,2));

Vstrike3 = [-R(2,3),R(1,3),0];
Vstrike3 = Vstrike3/norm(Vstrike3);
strike3 = atan2(Vstrike3(1),Vstrike3(2))*180/pi;
if isnan(strike3)
    strike3 = 0;
end
dip3 = acosd(R(3,3));

% calculates contribution of the first PTD
if DVx ~= 0
    [ue1,un1,uv1] = PTDdispSurf(X,Y,depth,strike1,dip1,DVx,nu);
else
    ue1 = zeros(size(X));
    un1 = zeros(size(X));
    uv1 = zeros(size(X));
end

% calculates contribution of the second PTD
if DVy ~= 0
    [ue2,un2,uv2] = PTDdispSurf(X,Y,depth,strike2,dip2,DVy,nu);
else
    ue2 = zeros(size(X));
    un2 = zeros(size(X));
    uv2 = zeros(size(X));
end

% calculates contribution of the third PTD
if DVz ~= 0
    [ue3,un3,uv3] = PTDdispSurf(X,Y,depth,strike3,dip3,DVz,nu);
else
    ue3 = zeros(size(X));
    un3 = zeros(size(X));
    uv3 = zeros(size(X));
end

% reshapes outputs to the original size of calculation points
sx = size(X);
ue = reshape(ue1 + ue2 + ue3,sx);
un = reshape(un1 + un2 + un3,sx);
uv = reshape(uv1 + uv2 + uv3,sx);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ue,un,uv]=PTDdispSurf(x,y,d,strike,dip,DV,nu)
% PTDdispSurf calculates surface displacements associated with a tensile 
% point dislocation (PTD) in an elastic half-space (Okada, 1985).

beta = strike - 90;
Rz = [cosd(beta) -sind(beta);sind(beta) cosd(beta)];
r_beta = Rz*[x y]';
x = r_beta(1,:)';
y = r_beta(2,:)';

r = (x.^2 + y.^2 + d.^2).^0.5;
q = y*sind(dip) - d*cosd(dip);

I1 = (1-2*nu)*y.*(1./r./(r+d).^2-x.^2.*(3*r+d)./r.^3./(r+d).^3);
I2 = (1-2*nu)*x.*(1./r./(r+d).^2-y.^2.*(3*r+d)./r.^3./(r+d).^3);
I3 = (1-2*nu)*x./r.^3 - I2;
I5 = (1-2*nu)*(1./r./(r+d)-x.^2.*(2*r+d)./r.^3./(r+d).^2);

% Note: For a PTD M0 = DV*mu!
ue = DV/2/pi*(3*x.*q.^2./r.^5-I3*sind(dip)^2);
un = DV/2/pi*(3*y.*q.^2./r.^5-I1*sind(dip)^2);
uv = DV/2/pi*(3*d.*q.^2./r.^5-I5*sind(dip)^2);

r_beta = Rz'*[ue un]';
ue = r_beta(1,:)';
un = r_beta(2,:)';

end
