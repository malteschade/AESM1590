function CMPnew=nmo_v(CMP, H_CMPgather, geo, c, smute)
% NMOedCMP = NMO_V(CMPgather, H_CMPgather, geo, c, smute) 
%
% This function applies NMO to a single CMP-gather, with linear
% interpolation, according to a constant velocity c.
%
% Output:     NMOedCMP - NMO-ed CMP-gather 
% Input:     CMPgather - CMP-gather
%          H_CMPgather - its header
%                  geo - geometry of the seismic data 
%                    c - constant velocity
%                smute - stretch-mute value (default is zero - no mute)
%
% See also NMO_VT, SORTDATA, SEMBLANCE

% Default the stretch-mute factor to zero
if ~exist('smute')
  smute = 0;
end

% Read the amount of time-samples and traces from the size of the datamatrix
[nt,nx]=size(CMP);

% Time sampling in [s]
dt = geo(1)/1000.;

% Initialise the CMP-gather after NMO
CMPnew = zeros(nt,nx);

if (smute == 0) 
  for ix=1:nx
      off    = H_CMPgather(2,ix);
      off2c2 = (off*off)/(c*c);
      for it=1:nt
        t0 = (it-1)*dt;
        t2 = t0*t0 + off2c2;
        tnmo = sqrt(t2) - t0;
        itnmo1    = floor(tnmo/dt);
        difft     = (tnmo-dt*itnmo1)/dt;
        if it+itnmo1+1 <= nt
            CMPnew(it,ix) = (1.-difft)*CMP(it+itnmo1,ix) + ...
                            difft*CMP(it+itnmo1+1,ix);
        end
        if it+itnmo1 == nt
            CMPnew(it,ix) = CMP(it+itnmo1,ix);
        end
      end
  end

else

  for ix=1:nx
    off    = H_CMPgather(2,ix);
    off2c2 = (off*off)/(c*c);
    for it=1:nt
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
                CMPnew(it,ix) = 0.0;
            else
                CMPnew(it,ix) = (1.-difft)*CMP(it+itnmo1,ix) + ...
                            difft*CMP(it+itnmo1+1,ix);
            end
        end
        if it+itnmo1 == nt
            CMPnew(it,ix) = CMP(it+itnmo1,ix);
        end
    end
  end

end 
