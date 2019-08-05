function varargout = plotpcdm(p,w,pview,varargin)
%PLOTPCDM Plot 3-D pCDM source
%   PLOTPCDM([X0,Y0,Z0,OX,OY,OZ,A,B],W) adds to current axe a 3-D
%   representation of pCDM source, with following parameters:
%		(X0,Y0,Z0) = center's coordinates
%		(OX,OY,OZ) = rotation angles (in degree) around X,Y and Z axis
%		A = horizontal over total volume variation ratio
%		B = vertical volume variation ratio
%		W = source maximum width
%
%	PLOTPCDM(...,'3d') plots the source in 3-D (default).
%	PLOTPCDM(...,'xy') plots a projection of the source in XY plan.
%	PLOTPCDM(...,'xz') plots a projection of the source in XZ plan.
%	PLOTPCDM(...,'zy') plots a projection of the source in ZY plan.
%
%
%	Authors: François Beauducel and Antoine Villié
%	Created: 2018-07-13, in Yogyakarta (Indonesia)
%	Updated: 2019-07-31



if nargin < 3
	pview = '3d';
end
opt = {'LineWidth',2,'EdgeColor','black','FaceColor','none'};

for n = 1:size(p,1)
	x0 = p(n,1);
	y0 = p(n,2);
	z0 = p(n,3);
	ox = p(n,4);
	oy = p(n,5);
	oz = p(n,6);
	a = p(n,7);
	b = p(n,8);

	% converts A and B aspect ratio to DVX/DVY/DVZ
	dvz = a;
	dvy = (1-dvz)*b;
	dvx = (1-dvz)*(1-b);

	dvmax = max([dvx,dvy,dvz]);
	dvz = w*dvz/dvmax;
	dvx = w*dvx/dvmax;
	dvy = w*dvy/dvmax;

	% unitary patch Z
	pz = [-1,-1, 1, 1;
	      -1, 1, 1,-1;
		   0, 0, 0, 0]/2;

	% unitary patch X
	px = [ 0, 0, 0, 0;
	      -1, 1, 1,-1;
		   1, 1,-1,-1]/2;

	% unitary patch Y
	py = [-1, 1, 1,-1;
		   0, 0, 0, 0;
	       1, 1,-1,-1]/2;

	% applies rotations
	Rx = [1 0 0;0 cosd(ox) sind(ox);0 -sind(ox) cosd(ox)];
	Ry = [cosd(oy) 0 -sind(oy);0 1 0;sind(oy) 0 cosd(oy)];
	Rz = [cosd(oz) sind(oz) 0;-sind(oz) cosd(oz) 0;0 0 1];
	R = Rz*Ry*Rx;

	px = R*(px*dvx);
	py = R*(py*dvy);
	pz = R*(pz*dvz);

	% plots result
	switch lower(pview)
		case 'xy'
			h(1) = patch(pz(1,:) + x0,pz(2,:) + y0,'k',opt{:},varargin{:});
			h(2) = patch(px(1,:) + x0,px(2,:) + y0,'k',opt{:},varargin{:});
			h(3) = patch(py(1,:) + x0,py(2,:) + y0,'k',opt{:},varargin{:});
		case 'xz'
			h(1) = patch(pz(1,:) + x0,pz(3,:) + z0,'k',opt{:},varargin{:});
			h(2) = patch(px(1,:) + x0,px(3,:) + z0,'k',opt{:},varargin{:});
			h(3) = patch(py(1,:) + x0,py(3,:) + z0,'k',opt{:},varargin{:});
		case 'zy'
			h(1) = patch(pz(3,:) + z0,pz(2,:) + y0,'k',opt{:},varargin{:});
			h(2) = patch(px(3,:) + z0,px(2,:) + y0,'k',opt{:},varargin{:});
			h(3) = patch(py(3,:) + z0,py(2,:) + y0,'k',opt{:},varargin{:});
			
		otherwise
			h(1) = patch(pz(1,:) + x0,pz(2,:) + y0, pz(3,:) + z0,'k',opt{:},varargin{:});
			h(2) = patch(px(1,:) + x0,px(2,:) + y0, px(3,:) + z0,'k',opt{:},varargin{:});
			h(3) = patch(py(1,:) + x0,py(2,:) + y0, py(3,:) + z0,'k',opt{:},varargin{:});
	end

	if nargout > 0
		varargout = h;
	end
	end
end
