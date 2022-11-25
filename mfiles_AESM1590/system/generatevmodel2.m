function vmodel = generatevmodel2(cmppicks, tvpicks, midpnts, geo)
% vmodel = GENERATEVMODEL2(cmppicks, tvpicks, midpnts, geo);
%
% This function builds a 2D velocity model of the subsurface, needed for
% the migration of the seismic data. The velocity model is built by linear
% interpolation from picked v(t)-logs. 
% After calculation the velocity model is plotted.
%
% Input: 
%   cmppicks - vector with CMP locations of time-velocity picks; size(ncmp)
%   tvpicks  - matrix with time-velocity picks; size(ncmp,npicks,2);
%              times in MILLIseconds, velocities in m/s
%   midpnts  - vector with all existing CMP-positions (from ANALYSEFOLD)
%   geo      - geometry of the seismic data
%
% Output: 
%     vmodel - velocity-model matrix (at all positions midpnts)
%
% Notes:
% - Make sure the sizes are consistent, so the number of cmp's (of vector 
%   cmppicks) corresponds to the number of rows of the matrix tvpicks
% - CMP numbers and times have to be in ascending order.
% - Matrix tvpicks is square so if you did not take the same number of 
%   time-velocity picks for the CMPs you took, you will have to copy 
%   velocities from your first timepick to earlier time(s) for those CMPs 
%   having fewer picks than the others; the velocity is assumed constant.
%
% Example of cmppicks vector and tvpicks matrix:
% cmppicks = [ 500; % 1st CMP position [m] for which you did velocity analysis
% 	          1000; % 2nd CMP position [m] for which you did velocity analysis
%             1500]
% tvpicks(:,:,1) = [  0     500  1000  1500; % picked times (ms) at 1st CMP
%                     0     400   900  1500; % picked times (ms) at 2nd CMP
%                     0     300   800  1500] % ...
% tvpicks(:,:,2) = [1500   1500  2000  2500; % picked velocities at 1st CMP
%                   1700   1700  2200  2700; % picked velocities at 2nd CMP
%                   1900   1900  2400  2900] % ...
%
% See also ANALYSEFOLD

% Read in the picked v(t)-logs
% velocitypicks;
t = tvpicks(:,:,1);
v = tvpicks(:,:,2);


% Calculating CMP sequence numbers [] from midpoint positions [m]
% Note that midpnts is an input argument from generatevmodel, the 
% function calling this script
cmpdist = midpnts(2)-midpnts(1);
cmp_initoffset = midpnts(1)/cmpdist;
for a=1:length(cmppicks);
    vcmp(a) = cmppicks(a)/cmpdist - (cmp_initoffset-1);
end

% Preparing the velocity and time matrices for generatevmodel.m
nrows=size(t,1);
ncols=size(t,2);
% Note that geo is an input argument from generatevmodel, the 
% function calling this script 
t_max=geo(1)*( geo(2) -1 );

% loop over rows, the amount of CMPs
for k=1:nrows
% Copy to zero timesample, otherwise copy zero-sample to first timesample
% in order to keep a square matrix
if t(k,1) > 0
    t_up(k,:) = [0      t(k,:)];
    v_up(k,:) = [v(k,1) v(k,:)];
else
    t_up(k,:) = [t(k,1) geo(1) t(k,2:ncols)];
    v_up(k,:) = [v(k,1) v(k,1) v(k,2:ncols)];
end

% Copy to t_max-timesample, otherwise copy t_max-sample to t=t_max-1 
% in order to keep a square matrix
% Beware that ncols of t_up has grown in size with one
if t_up(k,ncols+1) < t_max
    t_down(k,:) = [t_up(k,:)           t_max];
    v_down(k,:) = [v_up(k,:) v_up(k,ncols+1)];
else
    t_down(k,:) = [t_up(k,1:ncols) t_max-geo(1)    t_max          ]; 
    v_down(k,:) = [v_up(k,1:ncols) v_up(k,ncols+1) v_up(k,ncols+1)];
end
end

t=t_down;
v=v_down;


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
figure;
imagesc(midpnts,0:geo(1):geo(1)*geo(2)-1,vmodel);
xlabel('CMP position [m]');
ylabel('two-way time [ms]');
title('velocity model');
colorbar;
