function [seis, H_CSP, geo] = segyread(filename);
% [Data, H_CSP, geo] = SEGYREAD(filename);
%
% This function reads in a SEGY formatted file (with name and
% location given by 'filename'), into a matrix.
%
% Description of the output:
% Data  - seismic data in matrix format (traces in columns, 
%         timesamples in rows)
%
% H_CSP - shot-sorted header matrix in 5 rows 
%         Convention for the header-row numbering:
%         1 = the unique trace number
%         2 = trace offset with respect to shot position
%         3 = shot position corresponding to the trace
%         4 = receiver position corresponding to the trace
%         5 = cmp midpoint corresponding to the trace
%
% geo   - geometry of the seismic data. Its output is a row vector
%         of length 7, with the following elements:
%        (1) time sampling dt [ms]
%        (2) number of time samples nt in a trace 
%        (3) number of shots ns
%        (4) number of receivers nr per shot
%        (5) distance between two subsequent shots dxs [m]
%        (6) distance between two subsequent receivers dxr [m]
%        (7) absolute x-coordinate of the first shot x0 [m]
%
% Note that this function is only a wrapper for the CREWES functions
% readsegy and altreadsegy. For information on authors and terms of
% use please refer to those respective matlab-source files.
%

if nargout > 1
% call CREWES readsegy.m to do the hard work
%
% Note: this readsegy.m does not preserve amplitudes for our
% particular set of seismic data, but it can read headers.
% Original:
% [seis,t,line_name,ntr,nrec,xs,ys,xr,yr,offs,...
%		selevs,relevs,sdepths,cdps]=readsegy(filename);
% Adapted:
[seis,t,line_name,ntr,nrec,xs,ys,xr,yr,offs,...
		selevs,relevs,sdepths,cdps]=readsegy2(filename);
    
% Create custom headers used in this course, only if requested
    disp(' ')
    disp('  Bringing headers in standard format...')
    H_CSP=makehdr(t,line_name,ntr,nrec,xs,ys,xr,yr,offs,...
                       selevs,relevs,sdepths,cdps);
    geo=geometry(t,line_name,ntr,nrec,xs,ys,xr,yr,offs,...
                       selevs,relevs,sdepths,cdps);
    disp(' ')
else
    % this one does preserve amplitudes in our datasets
    % but does not always read headers correctly.
    [seis, sampint, textheader] = altreadsegy(filename, 'textheader','yes');
end