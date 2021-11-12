function s = pcdmdesc(a,b,tol)
%PCDMDESC pCDM source description
%	PCDMDESC(A,B) returns a string describing the point Compound Dislocation
%	Model source shape and orientation from the parameters:
%	   A: horizontal over total volume variation ratio A = dVZ/(dVX+dVY+dVZ)
%	   B: vertical volume variation ratio B = dVY/(dVX+dVY)
%
%	PCDMDESC(A,B,TOL) will use tolerancy TOL (TOL > 0 and < 1) to fit the
%	features instead of default 10% (TOL = 0.1).
%
%	Examples:
% 		A = 0, B = 0.5 : vertical pipe
% 		A = 0, B = 0 or 1 : vertical dyke
% 		A = 1/3, B = 0.5 : isotrop source
% 		A = 1, any B : horizontal sill
% 		A > 0, B > 0 : ellipsoid
%
%
%	Authors: F. Beauducel, A. Villié / WEBOBS
%	Created: 2019-02-21 in Yogyakarta (Indonesia)
%	Updated: 2020-04-17

if nargin < 2 || any([a,b] < 0 | [a,b] > 1)
	error('A and B must be scalars between 0 and 1');
end
if nargin < 3
	tol = 0.1;
end
if tol < 0 || tol > 0.5
	error('TOL must be scalar between 0 and 0.5');
end

s = 'ellipsoid';

if a < tol
	s = 'vert. ellipsoid';
	if b < tol
		s = 'vert. EW dyke';
	end
	if b > (1 - tol)
		s = 'vert. NS dyke';
	end
	if abs(b - 1/2) < tol
		s = 'vert. pipe';
	end
end

if abs(b - 1/2) < tol
	if abs(a - 1/3) < tol
		s = 'isotropic';
	end
	if a > (1/3 + tol)
		s = 'oblate ellipsoid';
	end
	if  a > tol && a < (1/3 - tol)
		s = 'prolate ellipsoid';
	end
end

if a > (1 - tol)
	s = 'horiz. sill';
end
