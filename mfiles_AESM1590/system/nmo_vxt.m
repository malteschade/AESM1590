function NMOedCMP=nmo_vxt(CMPgather, H_CMPgather, geo, c, smute)
% ==== system-script, do not execute manually ====
%
% NMOedCMP = NMO_VXT(CMPgather, H_CMPgather, geo, c, smute)
%
% Applies NMO to a single CMP gather with linear interpolation,
% according to a 1D velocity-time log c. In this case the VT-log
% is taken from the outcome of GENERATEVMODEL.
%
% This version is used by NMO_STACK. Use NMO_VT instead.
%
% Output:    NMOedCMP - NMO-ed CMP gather 
% Input:    CMPgather - CMP-gather
%         H_CMPgather - its header
%                   c - 1D velocity-time log at the CMP-position
%               smute - stretch-mute value (default is zero - no mute)

% Default the stretch-mute factor to zero
if ~exist('smute') 
  smute = 0;
end

nt   = geo(2);
dt   = geo(1)/1000.; % time in seconds
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
