function [gather_positions, gather_folds] =analysefold(H_SHT, sortkey)
% [positions, folds] = ANALYSEFOLD(H_CSP, sortkey);
%
% This function gives the positions of gathers, such as CMP's or 
% Common-Offsets gathers, as well as their folds. Furthermore, a
% crossplot is made.
%
% Clicking on a point displays the coordinates in the plot,
% <CTRL>-clicking makes them vanish again.
%
% Input:       H_CSP - shot-sorted data header
%            sortkey - sorting key, see below
% Output:  positions - vector with gather-positions
%              folds - vector with gather-folds
%
% Analysefold analyzes the fold of a dataset according to the 
% value of 'sortkey':
%    2 = Common Offset    (COF)
%    3 = Common Shot      (CSP)
%    4 = Common Receiver  (CRP)
%    5 = Common MidPoint  (CMP)
%
% See also GENERATEVMODEL

% Read the amount of time-samples and traces from the size of the datamatrix
[nt,ntr]=size(H_SHT);

% sort the header
H_SRT=sorthdr(H_SHT, sortkey);

% midpnt initialisation, 1st row contains midpnt distance, 2nd row the fold
out(1,1) = H_SRT(sortkey,1);
out(2,1) = 1;

% gather trace counter initialisation
% l is the amount of traces in a gather
l=1;

% Distance counter initialisation
% m is the amount of distances in the sorted dataset
m=1;

for k=2:ntr
	if H_SRT(sortkey,k) == out(1,m)
		l = l+1;
	else
		out(2,m)=l;
		m = m+1;
		out(1,m) = H_SRT(sortkey,k);
		l = 1;
	end
end

% Remove last superfluous column
out2=out(:,1:(m-1));

% Place in vectors for output
gather_positions = out2(1,:);
gather_folds = out2(2,:);

% Make a plot
figure;

plot(out2(1,:),out2(2,:),'x');
title('Amount of traces per gather');
xlabel('gather-distance [m]');
ylabel('fold');

datalabel;