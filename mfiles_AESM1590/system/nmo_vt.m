function NMOedCMP=nmo_vt(CMPgather, H_CMPgather, geo, smute)
% NMOedCMP = NMO_VT(CMPgather, H_CMPgather, geo, smute) 
%
% This function applies NMO to a single CMP-gather, with linear
% interpolation, according to a 1D velocity-time log: 
%
%    FIRST edit vlog.m, the input 1D velocity log 
%
% After calculation, three plots are made, showing the log, 
% and the CMP-gather before and after applying NMO.
%
%
% Output:     NMOedCMP - NMO-ed CMP-gather 
% Input:     CMPgather - CMP-gather
%          H_CMPgather - its header
%                  geo - geometry of the seismic data 
%                smute - stretch-mute value (default is zero - no mute)
%
% See also NMO_V, SORTDATA, SEMBLANCE

% Default the stretch-mute factor to zero
if ~exist('smute') 
  smute = 0;
end


% read in the used defined v-t log
c = vlog(geo);

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

