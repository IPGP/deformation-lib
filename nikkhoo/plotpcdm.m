function varargout = plotpcdm(p,w,varargin)
%PLOTPCDM Plot 3-D pCDM source
%   PLOTPCDM([X0,Y0,Z0,OX,OY,OZ,A,B],W) adds to current axe a 3-D
%   representation of pCDM source, with following parameters:
%		(X0,Y0,Z0) = center's coordinates
%		(OX,OY,OZ) = rotation angles (in degree) around X,Y and Z axis
%		A = horizontal over total volume variation ratio
%		B = vertical volume variation ratio
%		W = source maximum width
%
%	PLOTPCDM(...,'ellipsoid') plots a 3D ellipsoid instead of 3 rectangular
%	plans. Needs ELLIPSOID3 function.
%
%	PLOTPCDM(...,'3d') plots the source in 3-D (default).
%	PLOTPCDM(...,'xy') plots a projection of the source in XY plan.
%	PLOTPCDM(...,'xz') plots a projection of the source in XZ plan.
%	PLOTPCDM(...,'zy') plots a projection of the source in ZY plan.
%
%
%	Authors: François Beauducel and Antoine Villié
%	         inspired by the plotCDM function by Mehdi Nikkhoo (2016)

%	Created: 2018-07-13, in Yogyakarta (Indonesia)
%	Updated: 2020-06-08

pview = '3d';
pplan = true;
hflag = ishold;
if ~hflag
	hold on
end

if nargin > 2
	kk = cellfun(@ischar,varargin);
	
	% view type input argument
	k = ismember(varargin(kk),{'3d','xy','xz','zy'});
	if any(k)
		pview = varargin{kk(k)};
		varargin(kk(k)) = [];
	end
	
	% source type input argument
	k = strcmp(varargin(kk),'ellipsoid');
	if any(k)
		pplan = false;
		varargin(kk(k)) = [];
	end	
end

if isscalar(w)
	w = repmat(w,size(p,1),1);
end

%opt = {'LineWidth',2,'EdgeColor','black','FaceColor','none'};
opt = {'LineWidth',.1,'EdgeColor','black'};

% default dislocation colors by Nikkhoo
dcolx = [255 192 0]/255;
dcoly = [0 176  80]/255;
dcolz = [0 112 192]/255;

for n = 1:size(p,1)
	x0 = p(n,1);
	y0 = p(n,2);
	z0 = p(n,3);
	ox = p(n,4);
	oy = p(n,5);
	oz = p(n,6);
	a = p(n,7);
	b = p(n,8);

	% converts A and B aspect ratio to semi-axis DVX/DVY/DVZ with unitary
	% volume of the source
	dvx = 1/((1-a)*(1-b));
	dvy = 1/((1-a)*b);
	dvz = 1/a;

	dvmax = max([dvx,dvy,dvz]);
	dvz = w(n)*dvz/dvmax;
	dvx = w(n)*dvx/dvmax;
	dvy = w(n)*dvy/dvmax;

	if pplan
		% unitary patches (4 full squares)
		i4 = [0,-1,-1,0,0,0,-1,-1,0,1,1,0,0,0,1,1,0]/2;
		j4 = [0,0,-1,-1,0,1,1,0,0,0,1,1,0,-1,-1,0,0]/2;
		z4 = zeros(size(i4));

		pz = [i4*dvx;j4*dvy;z4];
		px = [z4;i4*dvy;j4*dvz];
		py = [i4*dvx;z4;j4*dvz];

		% applies rotations
		Rx = [1 0 0;0 cosd(ox) sind(ox);0 -sind(ox) cosd(ox)];
		Ry = [cosd(oy) 0 -sind(oy);0 1 0;sind(oy) 0 cosd(oy)];
		Rz = [cosd(oz) sind(oz) 0;-sind(oz) cosd(oz) 0;0 0 1];
		R = Rz*Ry*Rx;

		px = R*px;
		py = R*py;
		pz = R*pz;

		% plots result
		switch lower(pview)
			case 'xy'
				h(1) = patch(pz(1,:) + x0,pz(2,:) + y0,dcolz,opt{:},varargin{:});
				h(2) = patch(px(1,:) + x0,px(2,:) + y0,dcolx,opt{:},varargin{:});
				h(3) = patch(py(1,:) + x0,py(2,:) + y0,dcoly,opt{:},varargin{:});
			case 'xz'
				h(1) = patch(pz(1,:) + x0,pz(3,:) + z0,dcolz,opt{:},varargin{:});
				h(2) = patch(px(1,:) + x0,px(3,:) + z0,dcolx,opt{:},varargin{:});
				h(3) = patch(py(1,:) + x0,py(3,:) + z0,dcoly,opt{:},varargin{:});
			case 'zy'
				h(1) = patch(pz(3,:) + z0,pz(2,:) + y0,dcolz,opt{:},varargin{:});
				h(2) = patch(px(3,:) + z0,px(2,:) + y0,dcolx,opt{:},varargin{:});
				h(3) = patch(py(3,:) + z0,py(2,:) + y0,dcoly,opt{:},varargin{:});
			otherwise
				h(1) = patch(pz(1,:) + x0,pz(2,:) + y0, pz(3,:) + z0,dcolz,opt{:},varargin{:});
				h(2) = patch(px(1,:) + x0,px(2,:) + y0, px(3,:) + z0,dcolx,opt{:},varargin{:});
				h(3) = patch(py(1,:) + x0,py(2,:) + y0, py(3,:) + z0,dcoly,opt{:},varargin{:});
		end
		
	else
		[x,y,z] = ellipsoid3(dvx/2,dvy/2,dvz/2,ox,oy,oz,50);
		if nargout == 3
			varargout{1} = x;
			varargout{2} = y;
			varargout{3} = z;
		else
			h = surf(x + x0,y + y0,z + z0,varargin{:});
		end
	end
	
	if nargout == 1
		varargout{1} = h;
	end
	
	if ~hflag
		hold off
	end

end
