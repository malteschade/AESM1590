function out=plotgeometry(H, trmin, trmax, sortkey1, sortkey2)
% PLOTGEOMETRY(H, trmin, trmax, sortkey1, sortkey2)
%
% This function plots the geometry of dataset. It takes as input the 
% header matrix 'H' from the dataset and the first and last tracenumber
% 'trmin' and 'trmax' to be plotted. If 'trmin' and 'trmax' are not 
% specified, they default to the minimum and maximum trace present 
% in the seismic dataset.
%
% Additionally, before plotting, sorting of the header can take place 
% according to the value of 'sortkey1' as primary sort. Within 'sortkey1'
% the header is sorted by the value of 'sortkey2'. (sortkey1 must be 
% different from sortkey2). 
% 
% Valid values for sortkey are:
%    2 = Common Offset
%    3 = Common Shot (default)
%    4 = Common Receiver
%    5 = Common MidPoint
%
% See also SEGYREAD

%fill in some defaults if not specified by user
if ~exist('trmin')
  trmin=1;
end
if ~exist('trmax')
  trmax=max(size(H));
end

% fill in some sorting defaults if not specified
if ~exist('sortkey1')
    sortkey1 = 3;
end
if ~exist('sortkey2') 
  if (sortkey1 == 2) 
    sortkey2 = 5;
  end
  if (sortkey1 == 3)
    sortkey2 = 4;
  end
  if (sortkey1 == 4)
    sortkey2 = 3;
  end
  if (sortkey1 == 5)
    sortkey2 = 2;
  end
end

% Sort these headers
H_SRT=sorthdr(H, sortkey1, sortkey2);

% Start plotting
figure;

%trace_xs
subplot(2,2,1);
plot(trmin:trmax,H_SRT(3,trmin:trmax),'rx');
axis tight;
title('shot positions');
xlabel('trace number');
ylabel('shot position [m]');

%trace_cmp
subplot(2,2,2);
plot(trmin:trmax,H_SRT(5,trmin:trmax),'rx');
axis tight;
title('common midpoint positions (CMPs)');
xlabel('trace number');
ylabel('CMP [m]');

%trace_xr
subplot(2,2,3);
plot(trmin:trmax,H_SRT(4,trmin:trmax),'rx');
axis tight;
title('receiver positions');
xlabel('trace number');
ylabel('receiver position [m]');

%trace_off
subplot(2,2,4);
plot(trmin:trmax,H_SRT(2,trmin:trmax),'rx');
axis tight;
title('offsets ( = xs - xr )');
xlabel('trace number');
ylabel('offset [m]');
