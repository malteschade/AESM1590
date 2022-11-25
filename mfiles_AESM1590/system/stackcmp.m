function stackedtrace = stackcmp(gather)
% ==== system-script, do not execute manually ====
%
% stackedtrace = STACKCMP(gather)
%
% This function stacks one NMO-ed CMP-gather. Output is one stacked trace
% for the midpoint position belonging to the CMP-gather. Is used by the 
% function NMO_STACK, use STACKPLOT instead.
%
% See also NMO_STACK, STACKPLOT.

% count the number of traces in the gather
cmpsize=min(size(gather));

if cmpsize > 1
	stackedtrace = sum(gather,2);
else
	stackedtrace = gather;
end
