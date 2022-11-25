function out=plothdr(H,trmin,trmax)
% PLOTHDR(H, trmin, trmax) - plots the header data.
%
% This function takes as input the Header matrix 'H' and the first and last
% tracenumber 'trmin' and 'trmax' to be plotted. If 'trmin' and 'trmax' are
% not specified, they default to the minimum and maximum trace present in 
% the seismic dataset.
%
% Clicking on a point displays the coordinates in the plot,
% <CTRL>-clicking makes them vanish again.

%
% See also MAKEHDR.

%fill in some defaults if not specified by user
if ~exist('trmin')
  trmin=1;
end
if ~exist('trmax')
  trmax=max(size(H));
end

% Transpose matrix H for convenience
I = H';

% Start plotting
figure;

%trace_xs
al(1)=subplot(2,2,1);
plot(trmin:trmax,I(trmin:trmax,3),'x');
axis tight;
title('shot positions');
xlabel('trace number');
ylabel('shot position [m]');

%trace_cmp
al(2)=subplot(2,2,2);
plot(trmin:trmax,I(trmin:trmax,5),'x');
axis tight;
title('common midpoint positions (CMPs)');
xlabel('trace number');
ylabel('CMP [m]');

%trace_xr
al(3)=subplot(2,2,3);
plot(trmin:trmax,I(trmin:trmax,4),'x');
axis tight;
title('receiver positions');
xlabel('trace number');
ylabel('receiver position [m]');

%trace_off
al(4)=subplot(2,2,4);
plot(trmin:trmax,I(trmin:trmax,2),'x');
axis tight;
title('trace offset');
xlabel('trace number');
ylabel('offset [m]');

% Linking x-axes of every subplot while zooming a single one
%  works only for Matlab 7 and up, so first check version
v=version;
if str2num(v(1)) > 6;
    linkaxes(al,'x');
end

datalabel;