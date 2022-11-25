function semblance(CMPgather,H_CMPgather,geo, cmin,cmax,velstep)
% SEMBLANCE(CMPgather, H_CMPgather, geo, cmin, cmax, velstep) 
%
% This code generates a semblance panel for wave propagation velocities
% ranging from cmin to cmax, with a stepsize of velstep. As input it takes
% a CMP-gather at a certain midpoint together with its header.
%
% Input:    CMPgather   - CMP-gather
%           H_CMPgather - its header
%           geo         - geometry of the seismic data 
%           cmin        - minimum velocity  [m/s]
%           cmax        - maximum velocity  [m/s]
%           velstep     - velocity stepsize [m/s]
%
% See also SELECTCMP, NMO_VT

%constants
nt = geo(2);
dt = geo(1);
cmpsize=min(size(CMPgather));
counter=0;

disp(' ')
disp([' Velocity stepsize is ', num2str(velstep), ' m/s .'])
disp(' ')

% notation: NMO(i,j)

%call function nmo.m for making the NMO dataset with velocity ci
for ci=cmin:velstep:cmax;
disp([' processing for NMO velocity of ', num2str(ci), ' m/s ...'])
NMO=nmo_v(CMPgather,H_CMPgather,geo, ci);

counter = counter + 1;

for i=1:nt;
	sum = 0;
	sum2 = 0;
	for j=1:cmpsize;
		sum=sum+NMO(i,j);
		sum2=sum2+NMO(i,j)^2;
	end
	sqsum=sum^2;
	semtr(i,counter)=cmpsize^(-1)*(sqsum/sum2);
end

end

disp(' ')
plotimage(semtr,1:dt:(dt*nt),cmin:velstep:cmax);
title('Semblance plot')
xlabel('NMO velocity [m/s]')
ylabel('time [ms]')
