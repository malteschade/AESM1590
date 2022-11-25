function dmig = t2d(tmig,vmodel,geo, dz)
% DMIG = t2d(TMIG, vmodel, geo, dz)
%
% T2D converts a time-migrated section to depth, using a given velocity
% model, by means of 1D vertical time-to-depth conversion.
%
% Output:  DMIG - the migrated section in depth.
% Input:   TMIG - the time migrated section,
%        vmodel - the velocity model
%           geo - geometry of the seismic data
%            dz - output depth sampling interval
%
% See also KIRK_MIG2, GENERATEVMODEL.

% Get time sampling in [s]
dt=geo(1)/1000;

% Read the amount of time-samples and traces from the size of the
% time migrated datamatrix
[nt,cmp]=size(tmig);

% loop over migrated output positions
for j=1:cmp
    % Read in a trace with velocities    
    tvcol=vmodel(:,cmp);
    disp([' processing trace ', num2str(j), ' ...'])
    %generate the time-depth curve tz (see TIME2DEPTH).
    %note that we have vertical two-way traveltimes on the
    %time migrated section and stacking velocities in the vmodel
    tmp=0;
    for i=1:nt
        tmp=tmp+(dt/2)*tvcol(i);
        tz(i,1)=((i-1)*dt); % the two-way time sample
        tz(i,2)=tmp;        % the corresponding depth
    end
    dmig(:,j)=time2depth(tmig(:,j),0:dt:(nt-1)*dt,tz,dz);
end
 
