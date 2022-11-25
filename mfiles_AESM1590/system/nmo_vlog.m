function NMOedCMP=nmo_vlog(CMPgather, H_CMPgather, geo, t, v, smute)
% NMOedCMP = nmo_vlog(CMPgather, H_CMPgather, geo, t, v, smute) 
%
% This function applies NMO to a single CMP-gather, according to a 
% 1D velocity-time log
%
% Three plots are made, showing the log, and the CMP-gather before and 
% after applying NMO.
%
% Input:   CMPgather   - CMP-gather
%          H_CMPgather - its header
%          geo         - geometry of the seismic data 
%          t           - time of velocity picks (in seconds)
%          v           - velocitypicks (in m/s)
%          smute       - stretch-mute value (default is zero - no mute)
% Output:  NMOedCMP - NMO-ed CMP-gather 
%
% See also NMO_V, SORTDATA, SEMBLANCE

% Default the stretch-mute factor to zero
if ~exist('smute') 
  smute = 0;
end

% convert time in seconds to milliseconds
t = t*1000; 

if t(1) > 0  % do we need to fill up to t = 0 ms
    t = [0 t];
    v = [v(1) v];
end
if t(length(t)) < geo(1)*( geo(2) - 1 )  % do we need to fill down to t_max
    t = [t geo(1)*(geo(2)-1)];
    v = [v v(length(v))];
end

% Plotting time-veoocity log:
t2 = 0:geo(1):geo(1)*(geo(2)-1);
v2 = pwlint(t,v,t2);
figure;
plot(v2,t2,v,t,'r*');flipy;
xlabel('velocity [m/s]');ylabel('two-way traveltime [ms]');

c = v2';

nt   = geo(2);
dt   = geo(1)/1000.;
nx   = min(size(CMPgather));

NMOedCMP = zeros(nt,nx);

if (smute == 0)
  for ix=1:nx
    off    = H_CMPgather(2,ix);
    for it=1:nt
        off2c2 = (off*off)/(c(it)*c(it));
        t0 = (it-1)*dt;
        t2 = t0*t0 + off2c2;
        tnmo = sqrt(t2) - t0;
        itnmo1    = floor(tnmo/dt);
        difft     = (tnmo-dt*itnmo1)/dt;
        if it+itnmo1+1 <= nt
            NMOedCMP(it,ix) = (1.-difft)*CMPgather(it+itnmo1,ix) + ...
                            difft*CMPgather(it+itnmo1+1,ix);
        end
        if it+itnmo1 == nt
            NMOedCMP(it,ix) = CMPgather(it+itnmo1,ix);
        end
    end
  end
else
  for ix=1:nx
    off    = H_CMPgather(2,ix);
    for it=1:nt
        off2c2 = (off*off)/(c(it)*c(it));
        t0 = (it-1)*dt;
        t02 = t0*t0;
        t2 = t02 + off2c2;
        tnmo = sqrt(t2) - t0;
        if it == 1 
            dtnmo = 1000.;
        else
            dtnmo = abs(sqrt(1+off2c2/t02))-1.;
        end 
        itnmo1    = floor(tnmo/dt);
        difft     = (tnmo-dt*itnmo1)/dt;
        if it+itnmo1+1 <= nt
            if dtnmo >= smute
                NMOedCMP(it,ix) = 0.0;
            else
                NMOedCMP(it,ix) = (1.-difft)*CMPgather(it+itnmo1,ix) + ...
                            difft*CMPgather(it+itnmo1+1,ix);
            end
        end
        if it+itnmo1 == nt
            NMOedCMP(it,ix) = CMPgather(it+itnmo1,ix);
        end
    end
  end
end

figure; plotseis(CMPgather);
title('ORIGINAL CMP gather');

figure; plotseis(NMOedCMP);
title('NMO-ed CMP gather');

