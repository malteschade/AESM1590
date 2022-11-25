function out = geometry(t,line_name,ntr,nrec,xs,ys,xr,yr,offs,...
    selevs,relevs,sdepths,cdps)
% === system-script, do not execute manually ===
%
% GEOMETRY - defines the geometry of the seismic data.
%
% Output is a row vector of length 7, with the following elements:
% (1) time sampling dt [ms]
% (2) number of time samples nt in a trace 
% (3) number of shots ns
% (4) number of receivers nr per shot
% (5) distance between two subsequent shots dxs [m]
% (6) distance between two subsequent receivers dxr [m]
% (7) absolute x-coordinate of the first shot x0 [m]
%
% Notes: The input dataset is assumed to be shot-sorted.  Source- and
%        receiver-positions in the input SEGY-TraceHeaders are assumed
%        to be indicated in [mm]. x0 is set to zero.

% Note: Input arguments are taken from CREWES readsegy.m

% scalco (shorthand taken from Seismic Un*x) is the factor with which
% the coordinates have to be devided (in case of a neg. scalco). It is
% not evaluated by CREWES readsegy. Therefore we need to set it here
% (the chosen value is valid for the tripli.segy dataset supplied with
% this course).
scalco=-1000;              % i.e. positions in the headers are in [mm],
                           % multiply with 1000 to get [m]

% count number of traces in the dataset
trmax=ntr*nrec;

% reserve space
out=zeros(1,7);

% time sampling dt [ms]
out(1) = ( t(2)-t(1) )*1000;

% number of time samples nt in a trace 
out(2) = length(t);

% number of receivers nr per shot
% The code below assumes input data is shot-sorted
if trmax == 1
    % if dataset consists of single trace only
    out(4)=1;
else
    
i=0; a=0; b=0;
while a == b,
    i=i+1;
    a=xs(i);
    b=xs(i+1);
    if i+1 == trmax
        % if we have one shot, a will never be < b
        i=1;
        break
    end
end
out(4) = i ;

end

% number of shots ns
out(3) = trmax/out(4);

% distance between two subsequent shots dxs [m]
if ( out(3) == 1 )  % for handling data with only one shot
    out(5)= 0;
else
    out(5) = ( xs(out(4)+1)-xs(1) )/(-scalco);
end

% distance between two subsequent receivers dxr [m]
if ( out(4) == 1 ) % for handling data with only 1 receiver
    out(6)= 0;
else
    out(6) = ( xr(2)-xr(1) )/(-scalco);
end

% absolute x-coordinate of the first shot x0 [m]
out(7) = 0;
