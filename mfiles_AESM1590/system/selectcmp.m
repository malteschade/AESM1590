function [out,out2]=selectcmp(CMPsorted,H_CMP, midpnt)
% [CMPgather, H_CMPgather] = SELECTCMP(CMP_sorted_data, H_CMP, midpoint)
%
% This function selects plots a CMP gather according to its midpoint
% position. Midpoints can be found with the function ANALYSEFOLD.
%
% Input:  CMP_sorted_data - CMP-sorted dataset
%                   H_CMP - CMP-sorted data header
%                midpoint - midpoint of the CMP-gather you want to plot
% Output:       CMPgather - the CMP-gather
%             H_CMPgather - its header
%
% See also SORTDATA, ANALYSEFOLD

% Read the amount of time-samples and traces from the size of the datamatrix
[nt,ntr]=size(CMPsorted);

% CMP-gather trace counter initialisation:
% l is the number of traces in a CMP gather
l=1;

% Scan the CMP-sorted dataset for traces with the correct midpoint and put
% those traces in a cmpgather.
for k=1:ntr
	if H_CMP(5,k) == midpnt
		cmpgather(:,l)=CMPsorted(:,k);
		cmpgather_hdr(:,l)=H_CMP(:,k);
		l = l+1;
	end
end

out=cmpgather;
out2=cmpgather_hdr;


% Plot the requested CMP-gather
figure;
plotseis(cmpgather);
title(['CMP gather for midpoint at ', num2str(midpnt), ' m']);
