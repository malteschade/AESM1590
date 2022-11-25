function vmodel = generatevmodel(midpnts, geo)
% vmodel = GENERATEVMODEL(midpnts, geo);
%
% This function builds a 2D velocity model of the subsurface, needed for
% the migration of the seismic data. The velocity model is built by linear
% interpolation from picked v(t)-logs stored in the file velocitypicks.m: 
% 
%    FIRST edit velocitypicks.m to get proper output
%
% After calculation the velocity model is plotted.
%
% Input:   midpnts - the vector with CMP-positions (see ANALYSEFOLD)
%              geo - geometry of the seismic data 
% Output:   vmodel - the velocity model matrix
%
% See also ANALYSEFOLD

% Read in the picked v(t)-logs
velocitypicks;

% count CMP midpoints
m=max(size(midpnts));

% initialise the vmodel, the columns contain the velocities for each CMP midpnt
vxt=zeros(geo(2),m);

% copy first and last picks to the edges of the vmodel
vcmp_extended(1)=1;
v_extended(1,:)=v(1,:);
t_extended(1,:)=t(1,:);
for l=1:max(size(vcmp))
	vcmp_extended(l+1) = vcmp(l);
	v_extended(l+1,:) = v(l,:);
	t_extended(l+1,:) = t(l,:);
end
vcmp_extended(max(size(vcmp))+2)=m;
v_extended(max(size(vcmp))+2,:)=v(max(size(vcmp)),:);
t_extended(max(size(vcmp))+2,:)=t(max(size(vcmp)),:);

% time sample counter r=1;
for r=1:geo(2)
  %generate vlogs at the picked CMP positions k
  %vxt contains only the picked columns
  for k=1:max(size(vcmp_extended))
  	vxt(:,vcmp_extended(k)) = vlog_exd(v_extended(k,:),t_extended(k,:),geo);
	vold(r,k)=vxt(r,vcmp_extended(k));
  end
  % horizontal interpolation between the picked CMP positions
  x=vcmp_extended';
  x2=1:1:m;v2=pwlint(x,vold(r,:)',x2);
  vmodel(r,:)=v2;
end

% plot the velocity model
figure
imagesc(midpnts,0:geo(1):geo(1)*geo(2)-1,vmodel);
xlabel('CMP position [m]');
ylabel('two-way time [ms]');
title('velocity model');
colorbar;
