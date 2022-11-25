function [trc,seqno,itr,irec,dt,offset,sdepth,selev,relev,...
			xs,ys,xr,yr,cdp]=segyin2(fid)

% [trc,seqno,itr,irec,dt,offset,sdepth,selev,relev,...
%		xs,ys,xr,yr,cdp]=segyin(fid)
%
% Note: if end-of-file encountered in read, trc is returned as
%	nan and all other values are []
%
% fid ... valid file id
%
% trc ... trace
% seqno ... sequential trace number in the dataset
% itr ... trace number in current gather
% irec ... current gather number in line
% dt ... trace sample rate in seconds
% offset ... source to receiver offset
% sdepth ... source depth
% selev ... surface elevation at source
% relev ... surface elevation at receiver
% xs ... inline coordinate for source
% ys ... crossline coordinate for source
% xr ... inline coordinate for receiver
% yr ... crossline coordinate for receiver
% cdp ... cumulative cdp gather counter for line
% 
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.

% Modified: 1998-09-28, DSF
%    - when an EOF is encountered, all return variables are still
%      initialized.  This avoids the error message:
%    Warning: One or more output arguments not assigned during call to 'segyin'

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewes@geo.ucalgary.ca
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

% nsamp=length(trc);

%read the trace header
[stuff,count]=fread(fid,4,'int');
if(count<4)
	trc=nan;
        seqno=nan; itr=nan; irec=nan; dt=nan;
	offset=nan; sdepth=nan; selev=nan; relev=nan;
	xs=nan; ys=nan; xr=nan; yr=nan; cdp=nan;
	return;
end
seqno=stuff(1); irec=stuff(3); itr=stuff(4);
stuff=fread(fid,4,'char');
cdp=fread(fid,1,'int');
stuff=fread(fid,4,'char');
code=fread(fid,1,'short');
stuff=fread(fid,6,'char');
stuff=fread(fid,4,'int');
offset=stuff(1);relev=stuff(2);selev=stuff(3);sdepth=stuff(4);
stuff=fread(fid,16,'char');
stuff=fread(fid,2,'short');
stuff=fread(fid,4,'int');
xs=stuff(1);ys=stuff(2);xr=stuff(3);yr=stuff(4);
lenunit=fread(fid,1,'short');
stuff=fread(fid,24,'char');
stuff=fread(fid,2,'short');
nsamp=stuff(1);dt=stuff(2)/1000000.;
%Note SEGY stores dt in microseconds
stuff=fread(fid,122,'char');

%trace header complete

%read in trace
trc=fread(fid,nsamp,'float');


