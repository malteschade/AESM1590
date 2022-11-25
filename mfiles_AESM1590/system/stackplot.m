function stack = stackplot(gather, geo)
% stack = STACKPLOT(gather, geo)
%
% This function stacks the traces in a gather, and makes a plot
% of the stacked trace.
%
% Input:      geo - geometry of the seismic data 
%          gather - the gather you want to stack, usually
%                   this is an NMO-ed CMP-gather
% Output:   stack - stacked trace

gathersize=min(size(gather));
stack=sum(gather,2)/gathersize;

figure;
clf;
t       = (1:geo(1):geo(1)*geo(2));
[h,hva] = wtva(stack ,t);
flipy;
title('stacked trace')
xlabel('amplitude')
ylabel('time [ms]')
