function [ur,uz,B] = sun69(varargin)
%SUN69	Deformation from penny-shaped crack in elastic half-space.
%	[Ur,Uz] = SUN69(R,H,A,V) or SUN69(R,H,A,P,E,nu) computes radial and
%	vertical displacements Ur and Uz on the free surface, due to a 
%	horizontal circular fracture formed	in a semi-infinite elastic medium,
%	with following variables:
%		  R: radial distance of observation,
%		  H: depth of the center of the source from the surface,
%		  A: radius of the source with the hydrostatic pressure,
%		  V: volume of injected material,
%		  P: change of the hydrostatic pressure in the crack.
%		  E: Young's modulus,
%		 nu: Poisson's ratio (default is 0.25 for isotropic medium).
%
%	[Ur,Uz,B] = SUN69(...) returns also the maximum separation of fracture
%	(vertical displacement) B.
%
%	Equations from Sun [1969], with approximation H/A >> 1. If H/A > 2, error
%	is about 2 to 3%; if H/A > 5, solution is almost perfect.
%
%	Notes:
%		- Equations are all vectorized, so variables R,H,A,V or P can be 
%		  vectors or N-D matrix of the same size, while any of them can
%		  be scalars.
%		- Convention: Uz > 0 = UP, H is depth so in -Z direction.
%		- Units should be constistent, e.g.: R, H, A, Ur and Uz in m imply
%		  V in m3; optional E and P in Pa, nu dimensionless.
%
%	Example for a 3D plot of exagerated deformed surface:
%	  [x,y] = meshgrid(-3:.1:3);
%	  [th,rho] = cart2pol(x,y);
%	  [ur,uz] = sun69(rho,1,.5,1e6,10e9,0.25);
%	  [ux,uy] = pol2cart(th,ur);
%	  ps = 5e4;
%	  surf(x+ux*ps,y+uy*ps,uz*ps), axis equal, light
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%	Created: 2010
%	Updated: 2012-04-06
%
%	References:
%	 Sun, R. J. (1969). Theoretical size of hydraulically induced horizontal
%		fractures and corresponding surface uplift in an idealized medium,
%		J. Geophys. Res., 74, 5995-6011.

%	Copyright (c) 2012, François Beauducel, covered by BSD License.
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

error(nargchk(4,6,nargin))

for ii = 1:nargin
	if ~isnumeric(varargin{ii})
		error('All input arguments must be numeric.')
	end
end

% to check if input arguments have compatible sizes, constructs a complex
% vector of sizes, then uses UNIQUE on variables that are not scalar
sz = complex(cellfun('size',varargin,1),cellfun('size',varargin,2));
if length(unique(sz(find(sz~=complex(1,1))))) > 1
	error('All inputs must be scalar or matrix of the same size.')
end

h = varargin{2};
a = varargin{3};

har = min(h(:))/max(a(:));

if har < 2
	warning('Depth H must be greater than 2*A to obtain valid results.')
elseif har < 5
	warning('Depth H must be greater than 5*A. Possible 2-3% error on the solution.')
end

% cartesian to polar coordinates
r = varargin{1};

% B is the maximum separation of fracture given by equation (9) p. 5999
% and equation (10a) p. 5999 to specify volume variation
if nargin > 4
	if nargin == 5
		nu = varargin{6};
	else
		nu = 0.25;
	end
	% SUN69(R,H,A,P,E,nu)
	B = 8*(1 - nu.^2).*varargin{4}.*a./(pi.*varargin{5});
else
	% SUN69(R,H,A,V)
	B = 3*varargin{4}./(2*pi*a.^2);
end

% alternative: equation (14) p. 6001 in their complex form
%R1 = sqrt(r.^2 + (h - i*a).^2);
%R2 = sqrt(r.^2 + (h + i*a).^2);
%uz = i*B./(2*a) .* ((R1 - R2) + i*a.*h.*(1./R1 + 1./R2));

D = r.^2 + h.^2 - a.^2;
theta = atan2(2*a.*h,D);
k = ((D./a.^2).^2 + (2*h./a).^2).^.5;
srk = k.^.5;
akcos2 = a.*srk.*cos(theta/2);
aksin2 = a.*srk.*sin(theta/2);
hkcos2 = h.*srk.*cos(theta/2);
ak2cos = a.*k.*cos(theta);

% equation (16) p. 6001
uz = B.*(srk.*sin(theta/2) - (h./(a.*srk)).*cos(theta/2));

% equation (17) p. 6001
ur = (B.*r.*h./a) .* ((a + aksin2)./((h + akcos2).^2 + (a + aksin2).^2) ...
		- (hkcos2 - aksin2 + ak2cos)./((hkcos2 - aksin2 + ak2cos).^2 ...
		+ (akcos2 + h.*srk.*sin(theta/2) + a.*k.*sin(theta)).^2));

