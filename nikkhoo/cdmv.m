function varargout=cdmv(varargin)
%CDMV	Compound Dislocation Model
%	[uE,uN,uV,DV]=CDMV(X,Y,DEPTH,OMEGAX,OMEGAY,OMEGAZ,AX,AY,AZ,OPEN,NU)
%	calculates the surface displacements and potency associated with a
%	source composed of three mutually orthogonal rectangular dislocations
%	in a half-space (Nikkhoo et al., 2017).
%
%	This model is able to approximate displacements due to various shapes
%	of source like isotropic, sill, dyke, pipe, or any ellipsoid in any 
%	orientations in space.
%
%	Equations have been fully vectorized so all input parameters can be
%	vectors or N-D matrix of the same number of elements, and any of them
%	can be a scalar with some restrictions described below. Output 
%	parameters will have the same size as X.
%
%
%	--- Input parameters
%
%	X and Y:
%	   Horizontal coordinates (East, North) of calculation points relative
%	   to source located at (0,0). If X is a scalar, all other input 
%	   parameters must be also scalars. Use repmat(X,...) to convert X to
%	   a vector/matrix if necessary.
%
%	DEPTH:
%	   Depth of the source of the CDM centroid. The depth must be a positive 
%	   value and have the same unit as X and Y.
%
%	OMEGAX, OMEGAY and OMEGAZ:
%	   Clockwise rotation angles about X, Y and Z axes, respectively, that 
%	   specify the orientation of the CDM in space, in degrees. The three
%	   angles must have the same number of elements, i.e., mixing of scalar
%	   and vector/matrix is not allowed.
%
%	AX, AY and AZ:
%	   Semi-axes of the CDM along the X, Y and Z axes, respectively, before
%	   applying the rotations. Must have the same unit as X and Y.
%
%	OPEN:
%	   The opening (tensile component of the Burgers vector) of the RDs 
%	   that form the CDM. Must be the same as the unit of X and Y.
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
%	   same unit as OPEN and the CDM semi-axes in inputs.
%
%	DV:
%	   Potency of the CDM. DV has the unit of volume, i.e. the unit of 
%	   displacements, opening and CDM semi-axes to the power of 3.
%
% 
%	--- Glossary
%
%	 CDM: Compound Dislocation Model
%	  RD: Rectangular Dislocation
%	EFCS: Earth-Fixed Coordinate System
%	RDCS: Rectangular Dislocation Coordinate System
%	      The origin of the RDCS is the RD centroid. The axes of the RDCS 
%	      are aligned with the strike, dip and normal vectors of the RD,
%	      respectively.
%	ADCS: Angular Dislocation Coordinate System
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
%	Updated: 2019-09-03

%	Copyright (c) 2016 Mehdi Nikkhoo
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

% inputs arguments are converted to column vectors (linear indexing)
X  = varargin{1}(:);
Y  = varargin{2}(:);
D  = varargin{3}(:);
OX = varargin{4}(:);
OY = varargin{5}(:);
OZ = varargin{6}(:);
% AX, AY, AZ are converted to full axes
AX = 2*varargin{7}(:);
AY = 2*varargin{8}(:);
AZ = 2*varargin{9}(:);
OP = varargin{10}(:);
% NU is converted to 1-2*NU
if nargin < 11
	NU = 0.5; % defaut nu = 0.25
else
	NU = 1 - 2*varargin{11}(:);
end

% normalization of sizes to allow any mixing of scalar/matrix
if isscalar(OX)
    OX = repmat(OX,size(X));
end
if isscalar(OY)
    OY = repmat(OY,size(X));
end
if isscalar(OZ)
    OZ = repmat(OZ,size(X));
end
if isscalar(AX)
    AX = repmat(AX,size(X));
end
if isscalar(AY)
    AY = repmat(AY,size(X));
end
if isscalar(AZ)
    AZ = repmat(AZ,size(X));
end

% coefficients from the 3-D matrix of rotation
R11 =  cosd(OY).*cosd(OZ);
R12 = -cosd(OY).*sind(OZ);
R13 =  sind(OY);
R21 =  cosd(OZ).*sind(OX).*sind(OY) + cosd(OX).*sind(OZ);
R22 = -sind(OX).*sind(OY).*sind(OZ) + cosd(OX).*cosd(OZ);
R23 = -sind(OX).*cosd(OY);
R31 = -cosd(OX).*sind(OY).*cosd(OZ) + sind(OX).*sind(OZ);
R32 =  cosd(OX).*sind(OY).*sind(OZ) + sind(OX).*cosd(OZ);
R33 =  cosd(OX).*cosd(OY);

% coordinates for each RD summits
P1 = [AY.*R21/2 + AZ.*R31/2, AY.*R22/2 + AZ.*R32/2, AY.*R23/2 + AZ.*R33/2 - D];
P2 = P1 - repmat(AY,1,3).*[R21, R22, R23];
P3 = P2 - repmat(AZ,1,3).*[R31, R32, R33];
P4 = P1 - repmat(AZ,1,3).*[R31, R32, R33];

Q1 = [-AX.*R11/2 + AZ.*R31/2, -AX.*R12/2 + AZ.*R32/2, -AX.*R13/2 + AZ.*R33/2 - D];
Q2 = Q1 + repmat(AX,1,3).*[R11, R12, R13];
Q3 = Q2 - repmat(AZ,1,3).*[R31, R32, R33];
Q4 = Q1 - repmat(AZ,1,3).*[R31, R32, R33];

R1 = [AX.*R11/2 + AY.*R21/2, AX.*R12/2 + AY.*R22/2, AX.*R13/2 + AY.*R23/2 - D];
R2 = R1 - repmat(AX,1,3).*[R11, R12, R13];
R3 = R2 - repmat(AY,1,3).*[R21, R22, R23];
R4 = R1 - repmat(AY,1,3).*[R21, R22, R23];

[ue1,un1,uv1] = RDdispSurf(X,Y,P1,P2,P3,P4,OP,NU);
[ue2,un2,uv2] = RDdispSurf(X,Y,Q1,Q2,Q3,Q4,OP,NU);
[ue3,un3,uv3] = RDdispSurf(X,Y,R1,R2,R3,R4,OP,NU);

ue = ue1+ue2+ue3;
un = un1+un2+un3;
uv = uv1+uv2+uv3;

% special cases (one dimension is zero)
k = (AX==0 & AY~=0 & AZ~=0);
ue(k) = ue1(k);
un(k) = un1(k);
uv(k) = uv1(k);

k = (AX~=0 & AY==0 & AZ~=0);
ue(k) = ue2(k);
un(k) = un2(k);
uv(k) = uv2(k);

k = (AX~=0 & AY~=0 & AZ==0);
ue(k) = ue3(k);
un(k) = un3(k);
uv(k) = uv3(k);

% special cases (two or three dimensions are zero)
k = all([AX,AY,AZ]==0,2) | all([AX,AZ]==0,2) | all([AX,AZ]==0,2) | all([AY,AZ]==0,2);
ue(k) = 0;
un(k) = 0;
uv(k) = 0;

% half-space solution: The CDM must be under the free surface!
kair = any([P1(:,3) P2(:,3) P3(:,3) P4(:,3) ...
            Q1(:,3) Q2(:,3) Q3(:,3) Q4(:,3) ...
            R1(:,3) R2(:,3) R3(:,3) R4(:,3)] > 0,2);
ue(kair) = NaN;
un(kair) = NaN;
uv(kair) = NaN;

% reshapes outputs to the original size of X
varargout{1} = reshape(ue,sx);
varargout{2} = reshape(un,sx);
varargout{3} = reshape(uv,sx);

% Calculate the CDM total potency (AX, AY and AZ were converted to full axes)
if nargout > 3
	varargout{4} = reshape((AX.*AY + AX.*AZ + AY.*AZ).*OP,sx);
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ue,un,uv]=RDdispSurf(X,Y,P1,P2,P3,P4,OP,nu2)
% RDdispSurf calculates surface displacements associated with a rectangular
% dislocation in an elastic half-space.

Vnorm = cross(P2-P1,P4-P1,2);
Vnorm = Vnorm./repmat(sqrt(sum(Vnorm.^2,2)),1,3);
bX = OP.*Vnorm(:,1);
bY = OP.*Vnorm(:,2);
bZ = OP.*Vnorm(:,3);

[u1,v1,w1] = AngSetupFSC(X,Y,bX,bY,bZ,P1,P2,nu2); % Side P1P2
[u2,v2,w2] = AngSetupFSC(X,Y,bX,bY,bZ,P2,P3,nu2); % Side P2P3
[u3,v3,w3] = AngSetupFSC(X,Y,bX,bY,bZ,P3,P4,nu2); % Side P3P4
[u4,v4,w4] = AngSetupFSC(X,Y,bX,bY,bZ,P4,P1,nu2); % Side P4P1

ue = u1 + u2 + u3 + u4;
un = v1 + v2 + v3 + v4;
uv = w1 + w2 + w3 + w4;

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X1,X2,X3]=CoordTrans(x1,x2,x3,A1,A2,A3)
% CoordTrans transforms the coordinates of the vectors, from
% x1x2x3 coordinate system to X1X2X3 coordinate system. The transformation
% whose columns A1, A2 and A3 are the unit base vectors of the x1x2x3. The
% coordinates of e1,e2 and e3 in A must be given in X1X2X3. The transpose
% of A (i.e., A') will transform the coordinates from X1X2X3 into x1x2x3.

X1 = A1(:,1).*x1 + A1(:,2).*x2 + A1(:,3).*x3;
X2 = A2(:,1).*x1 + A2(:,2).*x2 + A2(:,3).*x3;
X3 = A3(:,1).*x1 + A3(:,2).*x2 + A3(:,3).*x3;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ue,un,uv]=AngSetupFSC(X,Y,bX,bY,bZ,PA,PB,nu2)
% AngSetupSurf calculates the displacements associated with an angular
% dislocation pair on each side of an RD in a half-space.

SideVec = PB-PA;
beta = acos(-SideVec(:,3)./sqrt(sum(SideVec.^2,2)));

ey1 = [SideVec(:,1:2),zeros(size(X))];
ey1 = ey1./repmat(sqrt(sum(ey1.^2,2)),1,3);
ey3 = repmat([0 0 -1],length(X),1);
%ey2 = cross(ey3,ey1,2);
ey2 = [ey1(:,2),-ey1(:,1),zeros(size(X))];

% Transform coordinates from EFCS to the first ADCS
[y1A,y2A,~] = CoordTrans(X-PA(:,1),Y-PA(:,2),-PA(:,3),ey1,ey2,ey3);
% Transform coordinates from EFCS to the second ADCS
[y1AB,y2AB,~] = CoordTrans(SideVec(:,1),SideVec(:,2),SideVec(:,3),ey1,ey2,ey3);
y1B = y1A-y1AB;
y2B = y2A-y2AB;

% Transform slip vector components from EFCS to ADCS
[b1,b2,b3] = CoordTrans(bX,bY,bZ,ey1,ey2,ey3);

[v1A,v2A,v3A] = AngDisDispSurf(y1A,y2A,beta,b1,b2,b3,nu2,-PA(:,3));
[v1B,v2B,v3B] = AngDisDispSurf(y1B,y2B,beta,b1,b2,b3,nu2,-PB(:,3));

% artefact-free for the calculation points near the free surface
I = (beta.*y1A)>=0;
if any(I)
    [v1A(I),v2A(I),v3A(I)] = AngDisDispSurf(y1A(I),y2A(I),beta(I)-pi,b1(I),b2(I),b3(I),nu2,-PA(I,3));
    [v1B(I),v2B(I),v3B(I)] = AngDisDispSurf(y1B(I),y2B(I),beta(I)-pi,b1(I),b2(I),b3(I),nu2,-PB(I,3));
end

% Calculate total displacements in ADCS
v1 = v1B - v1A;
v2 = v2B - v2A;
v3 = v3B - v3A;

% Transform total displacements from ADCS to EFCS
ue = ey1(:,1).*v1 + ey2(:,1).*v2;
un = ey1(:,2).*v1 + ey2(:,2).*v2;
uv = ey3(:,3).*v3;

k = abs(beta)<eps | abs(pi-beta)<eps;
ue(k) = 0;
un(k) = 0;
uv(k) = 0;

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v1,v2,v3] = AngDisDispSurf(y1,y2,beta,b1,b2,b3,nu2,a)
% AngDisDispSurf calculates the displacements associated with an angular
% dislocation in a half-space.

sinB = sin(beta);
cosB = cos(beta);
cotB = cot(beta);
z1 = y1.*cosB + a.*sinB;
z3 = y1.*sinB - a.*cosB;
r = sqrt(y1.^2 + y2.^2 + a.^2);

Fi = 2*atan2(y2,(r+a).*cot(beta/2) - y1); % The Burgers function

v1b1 = b1.*((1 - nu2.*cotB.^2).*Fi + y2./(r+a).*(nu2.*(cotB + y1/2./(r+a)) - y1./r) - y2.*(r.*sinB - y1).*cosB./r./(r-z3));
v2b1 = b1.*(nu2.*((.5+cotB.^2).*log(r+a)-cotB./sinB.*log(r-z3)) - 1./(r+a).*(nu2.*(y1.*cotB - a/2 - y2.^2./2./(r+a)) + y2.^2./r) + y2.^2.*cosB./r./(r-z3));
v3b1 = b1.*(nu2.*Fi.*cotB + y2./(r+a).*(1-nu2 + a./r) - y2.*cosB./(r-z3).*(cosB + a./r));

v1b2 = b2.*(-nu2.*((.5-cotB.^2).*log(r+a) + cotB.^2.*cosB.*log(r-z3))-1./(r+a).*(nu2.*(y1.*cotB + .5*a + y1.^2./2./(r+a))-y1.^2./r) + z1.*(r.*sinB-y1)./r./(r-z3));
v2b2 = b2.*((1 + nu2.*cotB.^2).*Fi-y2./(r+a).*(nu2.*(cotB + y1/2./(r+a))- y1./r) - y2.*z1./r./(r-z3));
v3b2 = b2.*(-nu2.*cotB.*(log(r+a) - cosB.*log(r-z3)) - y1./(r+a).*(1-nu2 + a./r) + z1./(r-z3).*(cosB + a./r));

v1b3 = b3.*(y2.*(r.*sinB - y1).*sinB./r./(r-z3));
v2b3 = b3.*(-y2.^2.*sinB./r./(r-z3));
v3b3 = b3.*(Fi + y2.*(r.*cosB + a).*sinB./r./(r-z3));

v1 = (v1b1 + v1b2 + v1b3)/2/pi;
v2 = (v2b1 + v2b2 + v2b3)/2/pi;
v3 = (v3b1 + v3b2 + v3b3)/2/pi;

end
