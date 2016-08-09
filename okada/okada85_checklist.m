function okada_checklist
% Checks consistency of OKADA85 results using the "Checklist for numerical calculations"
% from Table 2 [Okada, 1985] page 1149.
%
% References:
%	   Okada Y., Surface deformation due to shear and tensile faults in a
%	      half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.
%
% Author: Francois Beauducel <beauducel@ipgp.fr>
% Date: 2011-03-07


fprintf('\n=========================================\n');
fprintf('         Okada      Computed    Diff');
fprintf('\n=========================================\n');

% case 2
x = 2; y = 3; d = 4; dip = 70; L = 3; W = 2;
fprintf('------------- Case 2 -------------\n');
test = testok(1,x,y,d,dip,L,W,'strike');
test = testok(2,x,y,d,dip,L,W,'dip') + test;
test = testok(3,x,y,d,dip,L,W,'tensile') + test;

% case 3
x = 0; y = 0; d = 4; dip = 90; L = 3; W = 2;
fprintf('------------- Case 3 -------------\n');
test = testok(4,x,y,d,dip,L,W,'strike') + test;
test = testok(5,x,y,d,dip,L,W,'dip') + test;
test = testok(6,x,y,d,dip,L,W,'tensile') + test;

% case 4 (uses D=6, DIP=90, and RAKE=180 to simulate DIP=-90)
x = 0; y = 0; d = 6; dip = 90; L = 3; W = 2;
fprintf('------------- Case 4 -------------\n');
test = testok(7,x,y,d,dip,L,W,'strike',180) + test;
test = testok(8,x,y,d,dip,L,W,'dip') + test;
test = testok(9,x,y,d,dip,L,W,'tensile') + test;

fprintf('\n=========================================\n');
fprintf('Checklist results: %1.0f%% OK (%d/81 misfit).',100*test/81,81-test);
fprintf('\n=========================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = testok(j,x,y,d,dip,L,W,mode,rake)

if nargin < 9
	rake = 0;
end

ok = 0;
comp  = {'ux','uy','uz','dux/dx','dux/dy','duy/dx','duy/dy','duz/dx','duz/dy'};
varok = {'uE','uN','uZ',  '-uEE',  '-uEN',  '-uNE',  '-uNN',   'uZE',   'uZN'};	% NOTE: minus sign for stress convention: POSITIVE = COMPRESSION
check = [-8.689E-3,-4.298E-3,-2.747E-3,-1.220E-3,+2.470E-4,-8.191E-3,-5.814E-4,-5.175E-3,+2.945E-4;	% case 2 - strike
	-4.682E-3,-3.527E-2,-3.564E-2,-8.867E-3,-1.519E-4,+4.057E-3,-1.035E-2,+4.088E-3,+2.626E-3;	% case 2 - dip
	-2.660E-4,+1.056E-2,+3.214E-3,-5.655E-4,+1.993E-3,-1.066E-3,+1.230E-2,-3.730E-4,+1.040E-2;	% case 2 - tensile
	        0,+5.253E-3,        0,        0,-1.864E-2,-2.325E-3,        0,        0,+2.289E-2;	% case 3 - strike
	        0,        0,        0,        0,+2.748E-2,        0,        0,        0,-7.166E-2;	% case 3 - dip
	+1.223E-2,        0,-1.606E-2,-4.182E-3,        0,        0,-2.325E-3,-9.146E-3,        0;	% case 3 - tensile
	        0,-1.303E-3,        0,        0,+2.726E-3,+7.345E-4,        0,        0,-4.422E-3;	% case 4 - strike
	        0,        0,        0,        0,+5.157E-3,        0,        0,        0,-1.901E-2;	% case 4 - dip
	+3.507E-3,        0,-7.740E-3,-1.770E-3,        0,        0,-7.345E-4,-1.843E-3,        0];	% case 4 - tensile

switch mode
case 'strike'
	slip = 1;
	u3 = 0;
case 'dip'
	rake = 90;
	slip = 1;
	u3 = 0;
case 'tensile'
	slip = 0;
	u3 = 1;
end
[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(x-L/2,y-cosd(dip)*W/2,d-sind(dip)*W/2,90,dip,L,W,rake,slip,u3);

fprintf('--- Mode %s\n',mode);
for i = 1:length(comp)
	v = okdiff(check(j,i),eval(varok{i}));
	fprintf('%7s  %+1.3e  %+1.3e  %g  ',comp{i},v);
	if v(3)==0
		fprintf('ok\n');
		ok = ok + 1;
	else
		fprintf('* misfit *\n');
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = okdiff(o,x)

d = roundsd(x,4)-o;
if abs(d) < eps
	d = 0;
end
if abs(x) < eps
	x = 0;
end
v = [o,x,d];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=roundsd(x,n)

og = 10.^(floor(log10(abs(x)) - n + 1));
y = round(x./og).*og;
y(x==0) = 0;