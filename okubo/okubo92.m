function varargout=okubo92(varargin)
%OKUBO92 Surface gravity change due to a finite rectangular source.
%	[dG,dH] = OKUBO92(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,RHO,RHOP)
%	computes total gravity change at the free surface of an elastic
%	half-space, due to a dislocation defined by RAKE, SLIP, and OPEN on a 
%	rectangular fault defined by orientation STRIKE and DIP, and size LENGTH
%	and WIDTH. The fault centroid is located (0,0,-DEPTH).
%
%	   E,N    : vector or matrix of coordinates of observation points in a 
%	            geographic referential (East,North,Up) relative to fault 
%	            centroid (units are described below)
%	   DEPTH  : depth of the fault centroid (DEPTH > 0)
%	   STRIKE : fault trace direction (0 to 360° relative to North), defined
%	            so that the fault dips to the right side of the trace
%	   DIP    : angle between the fault and a horizontal plane (0 to 90°)
%	   LENGTH : fault length in the STRIKE direction (LENGTH > 0)
%	   WIDTH  : fault width in the DIP direction (WIDTH > 0)
%	   RAKE   : direction the hanging wall moves during rupture, measured 
%	            relative to the fault STRIKE (-180 to 180°).
%	   SLIP   : dislocation in RAKE direction (in m)
%	   OPEN   : dislocation in tensile component (in m)
%	   RHO    : density of the medium (in kg/m^3)
%	   RHOP   : optional density of the cavity-filling matter (in kg/m^3); 
%	            takes RHO value if not specified.
%
%	returns the following variables (same matrix size as E and N):
%	   dG : total gravity change (in m/s^2)
%	   dH : elevation changes as given by Okada [1985] model (in m)
%
%	Length unit consistency: E, N, DEPTH, LENGTH, and WIDTH must have the  
%	same unit (e.g. km) which can be different from that of SLIP and OPEN
%	(e.g. m to to produce dG in m/s^2 and dH in m).
%
%	[...] = OKUBO92(...,RHOP,BETA,NU) specifies gravity vertical gradient  
%	BETA (default is 0.309e-5 s^-2) and Poisson's ratio NU (default is 0.25 
%	for an isotropic medium).
%
%	Formulas and notations from Okubo [1992] solution excepted for the fault 
%	geometry parameters after Aki & Richards [1980], e.g.:
%	      DIP=90, RAKE=0   : left lateral (senestral) strike slip
%	      DIP=90, RAKE=180 : right lateral (dextral) strike slip
%	      DIP=70, RAKE=90  : reverse fault
%	      DIP=70, RAKE=-90 : normal fault
%
%	Equations are all vectorized excepted for argument DIP which must be
%	a scalar; all other arguments can be scalar or matrix of the same size.
%
%	Example:
%
%	   [E,N] = meshgrid(linspace(-10,10,50));
%	   dG = okubo92(E,N,6,90,90,10,10,0,5,0,2670);
%	   figure, pcolor(E,N,1e8*dG), shading interp
%	   hold on, [c,h]=contour(E,N,1e8*dG,-50:10:50,'k'); clabel(c,h), hold off
%	   colorbar, polarmap
%
%	reproduces the figure 4a in Okubo [1992] paper, considering a 10x10km 
%	fault at 6km depth, oriented N90°-strike, 90°-dip, and 5m senestral 
%	dislocation. Gravity changes are computed on a regular grid from -10 to  
%	10km, and are plotted as a surface in 1e-8 m/s^2 unit. POLARMAP is an 
%	author's colormap.
%
%	See also OKADA85 author's function to compute complete displacements,
%	tilt and strain at free surface.
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%	   Institut de Physique du Globe de Paris
%	Created: 2012-06-11
%	Updated: 2013-07-26
%
%	References:
%	   Aki K., and P. G. Richards, Quantitative seismology, Freemann & Co,
%	      New York, 1980.
%	   Okubo S., Gravity and Potential Changes due to Shear and Tensile Faults
%	      in a Half-Space, J. Geophys. Res., 97:B5, 7137-7144, 1992.
%
%	Acknowledgments: Bernard F. Whiting, Eric Chassande-Mottin, Martin
%	Fuchs

%	Copyright (c) 2013, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

if nargin < 11 || nargin > 14
	error('Wrong number of input arguments.')
end

for i = 1:nargin
	if ~isnumeric(varargin{i})
		error('Input arguments E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,RHO,RHOP must be all numeric.')
	end
end
if numel(varargin{5}) ~= 1
	error('DIP argument must be scalar.')
end

% Newton's gravitational constant
G = 6.67384e-11;	% (m^3/kg/s^2)

% Default values for optional input arguments
nu = 0.25;	% isotropic Poisson's ratio
beta = 0.309e-5;	% free-air gravity gradient (m.s^-2/m)

% Assigns input arguments
e = varargin{1};
n = varargin{2};
depth = varargin{3};
strike = varargin{4};
dip = varargin{5};	% 'delta' in Okubo's equations
L = varargin{6};
W = varargin{7};
rake = varargin{8};
slip = varargin{9};
U3 = varargin{10};
rho = varargin{11};

if nargin > 11
	rhop = varargin{12};
else
	rhop = rho;
end

if nargin > 13
	beta = varargin{13};
	nu = varargin{14};
end

% Defines dislocation in the fault plane system
U1 = cosd(rake).*slip;
U2 = sind(rake).*slip;

% Converts fault coordinates (E,N,DEPTH) relative to centroid
% into Okada's reference system (X,Y,D)
d = depth + sind(dip).*W/2;	% d is fault's base edge as Okada's convention
ec = e + cosd(strike)*cosd(dip).*W/2;
nc = n - sind(strike)*cosd(dip).*W/2;
x = cosd(strike).*nc + sind(strike).*ec + L/2;
y = sind(strike).*nc - cosd(strike).*ec + cosd(dip).*W;

% Variable substitution (independent from xi and eta)
p = y.*cosd(dip) + d.*sind(dip);
q = y.*sind(dip) - d.*cosd(dip);

% Elevation changes dH (must be computed first) [equation (57) p. 7139]
dH = 1/(2*pi) * (U1*chinnery(@Sh,x,p,L,W,q,dip,nu) ... % strike-slip
	+ U2 * chinnery(@Dh,x,p,L,W,q,dip,nu) ... % dip-slip
	+ U3 * chinnery(@Th,x,p,L,W,q,dip,nu) ... % tensile fault
	);

% Total gravity changes dG [equation (49) p. 7139]
dG = rho*G * (U1*chinnery(@Sg,x,p,L,W,q,dip,nu) ... % strike-slip
	+ U2 * chinnery(@Dg,x,p,L,W,q,dip,nu) ... % dip-slip
	+ U3 * chinnery(@Tg,x,p,L,W,q,dip,nu)) ... % tensile fault
     + (rhop-rho)*G*U3*chinnery(@Cg,x,p,L,W,q,dip,nu) ... % filling cavity
     - beta*dH; % free-air effect

% Assigns output arguments
switch nargout
	case 1
		varargout = {dG};
	case 2
		varargout = {dG,dH};
	otherwise
		disp('Unvalid number of output arguments.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =================================================================
% --- Chinnery's notation [equation (23) p. 7138]
function u=chinnery(f,x,p,L,W,q,dip,nu)
u = feval(f,x,p,q,dip,nu) ...
	- feval(f,x,p-W,q,dip,nu) ...
	- feval(f,x-L,p,q,dip,nu) ...
	+ feval(f,x-L,p-W,q,dip,nu);


% =================================================================
% Displacement subfunctions (Okada model regiven by Okubo)

% --- strike-slip displacement subfunction [equation (58) p. 7139]
function u=Sh(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sind(dip) - q*cosd(dip);
u = - db.*q./(R.*(R + eta)) ...
	- q.*sind(dip)./(R + eta) ...
	- I4(db,eta,q,dip,nu,R).*sind(dip);

% --- dip-slip displacement subfunction [equation (59) p. 7139]
function u=Dh(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sind(dip) - q*cosd(dip);
u = - db.*q./(R.*(R + xi)) ...
	- sind(dip).*atan(xi.*eta./(q.*R)) ...
	+ I5(xi,eta,q,dip,nu,R,db).*sind(dip).*cosd(dip);

% --- tensile fault displacement subfunction [equation (60) p. 7139]
function u=Th(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sind(dip) - q*cosd(dip);
u = (eta*cosd(dip) + q*sind(dip)).*q./(R.*(R + xi)) ...
	+ cosd(dip).*(xi.*q./(R.*(R + eta)) ...
	- atan(xi.*eta./(q.*R))) ...
	- I5(xi,eta,q,dip,nu,R,db).*sind(dip).^2;

% --- I4 displacement subfunction [equations (61) and (63) p. 7139-7140]
function I=I4(db,eta,q,dip,nu,R)
if cosd(dip) == 0
	I = -(1 - 2*nu) * q./(R + db);
else
	I = (1 - 2*nu) * 1/cosd(dip) * (log(R + db) - sind(dip)*log(R + eta));
end

% --- I5 displacement subfunction [equations (62) and (64) p. 7140]
function I=I5(xi,eta,q,dip,nu,R,db)
%X = sqrt(xi.^2 + q.^2);
if cosd(dip) == 0
	I = -(1 - 2*nu) * xi.*sind(dip)./(R + db);
else
	I = (1 - 2*nu) * 2./cosd(dip) ...
		.* atan((-q.*cosd(dip) + (1 + sind(dip)).*(R + eta))./(xi.*cosd(dip)));
%		.* atan((eta.*(X + q.*cosd(dip)) + X.*(R + X).*sind(dip))./(xi.*(R + X).*cosd(dip)));
end


% =================================================================
% Gravity subfunctions

% --- strike-slip gravity subfunction [equation (52) p. 7139]
function u=Sg(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = - q*sind(dip)./R + (q.^2 * cosd(dip))./(R.*(R + eta));

% --- dip-slip gravity subfunction [equation (53) p. 7139]
function u=Dg(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sind(dip) - q*cosd(dip);
u = 2*I2(xi,eta,q,R)*sind(dip) - q.*db./(R.*(R + xi));

% --- tensile fault gravity subfunction [equation (54) p. 7139]
function u=Tg(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta*cosd(dip) + q*sind(dip);
u = 2*I2(xi,eta,q,R)*cosd(dip) + q.*yb./(R.*(R+xi)) + q.*xi*cosd(dip)./(R.*(R+eta));

% --- cavity-filling gravity subfunction [equation (55) p. 7139]
function u=Cg(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = 2*I2(xi,eta,q,R)*cosd(dip) - sind(dip)*log(R + xi);

% --- I2 subfunction [equation (32) p. 7139]
function I=I2(xi,eta,q,R)
I = atan((R + xi + eta)./q);

