function zosection=nmo_stack(cmpsorted,cmpsorted_hdr, midpnts, ...
    folds, geo, vmodel, smute)
% zosection = NMO_STACK(cmpsorted, cmpsorted_hdr, midpnts, ...
%                                         folds, geo, vmodel, smute);
%
% This function generates a stacked zero-offset section from a CMP-sorted
% dataset. First, NMO correction is performed on each CMP-gather, using the
% velocity model of the subsurface. Subsequently, each NMO'ed CMP-gather is
% stacked to a obtain zero-offset traces on the distances corresponding to 
% the midpoints of each CMP-gather.
%
% Input:      cmpsorted - CMP-sorted dataset
%         cmpsorted_hdr - its headers
% 	            midpnts - vector with CMP-gather positions (see ANALYSEFOLD)
% 	              folds - vector with CMP-gather folds (see ANALYSEFOLD)
%                   geo - geometry of the seismic data 
%                vmodel - velocity model matrix
%                 smute - stretch-mute factor (default is 0 - no mute)
% Output:     zosection - zero-offset stacked seismic section
%
% See also SORTDATA, ANALYSEFOLD, GENERATEVMODEL

% Default the stretch-mute factor to zero
if ~exist('smute') 
  smute = 0;
end

% Read the amount of time-samples and traces from the size of the datamatrix
[nt,ntr]=size(cmpsorted);

% Amount of CMP-gathers equals the length of CMP midpoint-vector
cmpnr=length(midpnts);

% Initialise tracenr in cmpsorted dataset
tracenr=1;

% For each CMP-gather, do ...
for l=1:cmpnr
    % CMP midpoint in [m] (just for display), and associated fold
    midpoint=midpnts(l);
    fold=folds(l);
    disp([' processing CMP at ', ...
            num2str(midpoint), ' m, with fold ', num2str(fold), ' ...'])
    % positioning in the cmpsorted dataset
    gather=cmpsorted(:,tracenr:tracenr+(fold-1) );
    gather_hdr=cmpsorted_hdr(:,tracenr:tracenr+(fold-1) );
    % NMO and stack the selected CMP-gather
        nmoed=nmo_vxt(gather,gather_hdr, geo, vmodel(:,l),smute);
		zotrace=stackcmp(nmoed);
		zosection(:,l)=zotrace;
    % go to traceposition of next CMP in cmpsorted dataset
    tracenr=tracenr+fold;
end
