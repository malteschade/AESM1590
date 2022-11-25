% VELOCITYPICKS - script containing velocity picks for selected CMPs.
%
% Each CMP must contain the same number of picks, five will generally do.
% For times smaller than the first time-pick or larger than the last 
% time-pick, the velocity is assumed constant.
%
% cmppos - CMP position [m]
% v      - NMO velocity [m/s]
% t      - time of picked NMO velocity [ms]


% --------- begin editing here ---------
%
% Given numbers are just examples.
% Times and CMP numbers have to be in ascending order.
% Matrices t and v have to be square: this means that if you did not take the 
%  same amount of velocity-time picks for all CMPs, you will have to copy 
%  velocities from your first timepick to earlier time(s), for those CMPs which
%  have fewer picks than the others (hereby assuming the velocity in the 
%  overburden stays constant).
%
cmppos = [ 800;  % the 1st CMP position [m] for which you did velocity analysis
	      1152;  % the 2nd CMP position [m] for which you did velocity analysis
          1504;  % etc ...
	      1848;
          2200]
t=[0 180 350 370 590 800 1400;   % time picks [ms] for 1st CMP
   0 180 350 370 590 800 1400;   % time picks [ms] for 2nd CMP
   0 180 350 370 590 800 1400; % etc ...
   0 370 500 580 800 1300 1400;
   0 370 500 580 800 1300 1400]
v=[1475 1475 1600 1800 1800 2050 2050; % velocities [m/s] at times picked for 1st CMP
   1475 1475 1600 1800 1800 2050 2050; % velocities [m/s] at times picked for 2nd CMP
   1600 1600 1600 2050 2050 2050 2050; % etc ...
   1450 1450 1800 1800 2150 2150 2150;
   1450 1450 1800 1800 2150 2150 2150]
%
% --------- end editing here -----------


% Calculating CMP sequence numbers [] from midpoint positions [m]
% Note that midpnts is an input argument from generatevmodel, the 
% function calling this script
cmpdist = midpnts(2)-midpnts(1);
cmp_initoffset = midpnts(1)/cmpdist;
for a=1:length(cmppos);
    vcmp(a) = cmppos(a)/cmpdist - (cmp_initoffset-1);
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
