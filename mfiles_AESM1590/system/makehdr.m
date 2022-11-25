function H_CSP=makehdr(t,line_name,ntr,nrec,xs,ys,xr,yr,offs,...
                       selevs,relevs,sdepths,cdps)
% === system script, do not execute manually ===
%
% MAKEHDR - reformat the headers from the raw 
%           (i.e. shot-sorted) seismic dataset.
%
% Output is H_CSP, the custom 5-row header matrix used in this course.
% Logically, the header has as much columns as the dataset has traces,
% so that each header-column corresponds to the corresponding trace 
% (column) in the seismic data-matrix.
%
% Convention for contents of each header-row:
%    row 1 = the unique trace numbers
%    row 2 = trace offsets with respect to shot position
%    row 3 = shot positions corresponding to each trace
%    row 4 = receiver positions corresponding to each trace
%    row 5 = cmp midpoints corresponding to each trace
%
% See also PLOTHDR.

% Note: Input arguments are taken from CREWES readsegy.m

% scalco (shorthand taken from Seismic Un*x) is the factor with which
% the coordinates have to be devided (in case of a neg. scalco) to arrive
% at meters. It is not evaluated by CREWES readsegy. Therefore we need
% to set it here (the chosen value is valid for the tripli.segy dataset
% supplied with this course).
scalco=-1000;              % i.e. positions in the headers are in [mm],
                           % devide by 1000 to get [m]

% Generate 5 headerrows:
%
% 1. the unique trace number
trace_nr=1:ntr*nrec;

% 2. trace offset with respect to shot position
trace_off=offs;

% 3. shot position at which this trace recording was made
trace_xs=xs/(-scalco);

% 4. receiver position of this trace recording
trace_xr=xr/(-scalco);

% 5. cmp midpoint corresponding to the trace
trace_cmp=cdps/(-scalco);

% Build the header matrix
H_CSP=zeros(5,ntr*nrec);
H_CSP(1,:)=trace_nr;
H_CSP(2,:)=trace_off;
H_CSP(3,:)=trace_xs;
H_CSP(4,:)=trace_xr;
H_CSP(5,:)=trace_cmp;
