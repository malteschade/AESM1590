function [Data,H_CSP,geo]=segyread(filename,endian,IBM)
% [Data, H_CSP, geo] = SEGYREAD(filename);
%
% This function reads in a SEGY formatted file (with name and
% location given by 'filename'), into a matrix.
%
% Description of the output:
% Data  - seismic data in matrix format (traces in columns, 
%         timesamples in rows)
%
% H_CSP - shot-sorted header matrix in 5 rows 
%         Convention for the header-row numbering:
%         1 = the unique trace number
%         2 = trace offset with respect to shot position
%         3 = shot position corresponding to the trace
%         4 = receiver position corresponding to the trace
%         5 = cmp midpoint corresponding to the trace
%
% geo   - geometry of the seismic data. Its output is a row vector
%         of length 7, with the following elements:
%        (1) time sampling dt [ms]
%        (2) number of time samples nt in a trace 
%        (3) number of shots ns
%        (4) number of receivers nr per shot
%        (5) distance between two subsequent shots dxs [m]
%        (6) distance between two subsequent receivers dxr [m]
%        (7) absolute x-coordinate of the first shot x0 [m]
%
% [Data, H_CSP, geo] = SEGYREAD(filename, endian) allows file endianess to
% be specified. According to SEGY Rev. 1 specifications all SEGY files
% should be IBM big endian. Specifying 'L' will indicate that the file is
% saved as 'ieee-l' small endian, which is common for desktop PCs. 'B' is
% default IBM big endian encoding of the SEGY file.
%
% NOTE: By using this software, you are agreeing to the terms detailed in 
% this software's Matlab source file. Use completely at your own risk.

% 2015: Speed greatly improved by Max Holicki

%% Check Inputs & Outputs
% endian
if nargin==1;
    endian='B';
    IBM=true;
elseif nargin==2;
    IBM=true;
elseif nargin==0||nargin>3;
    error('segyread: Too many or too few input arguments.');
elseif ~ischar(endian)||numel(endian)~=1;
    error('segyread: Incorrect endian specification.');
end;
% filename
if ~ischar(filename);
    error('segyread: filename is not a string.');
end;
% nargout
if nargout>3;
    error('segyread: too many outputs.');
end;
%% Get Endianess & IBM Float Conversion
IBM=(endian=='B'&&IBM);
[~,~,E]=computer;
E=(E=='L'&&endian=='B')||(E=='B'&&endian=='L');
%% Read File to Buffer
fid=fopen(filename);
if fid<0;
    error('segyread: filename is not a valid filename or path.');
end;
buff=fread(fid,'*uint8');
fclose(fid);
%% Extract Important SEGY Headers
dt=typecast(buff(3217:3218),'uint16');
nt=typecast(buff(3221:3222),'uint16');
% Swap Endianess
if E;
    dt=swapbytes(dt);
    nt=swapbytes(nt);
end;
%% Allocate Traces
NT=uint64(nt);
NTRC=(numel(buff)-3840)/(240+NT*4);
i=3841+bsxfun(@plus,(0:NTRC-1).*(240+NT*4),uint64(0:nt*4-1)');
Data=reshape(typecast(buff(i(:)),'SINGLE'),NT,NTRC);
% Swap Endianess
if E;
    Data=swapbytes(Data);
end;
% Change Float Type
if IBM;
    %tmp=typecast(Data(:),'uint32');
    %for i=1:numel(Data);Data(i)=ibm2ieee(tmp(i));end;
    Data=reshape(ibm2ieee(typecast(Data(:),'uint32')),size(Data));
end;
%% Allocate Trace Headers
if nargout>1;
    i=3600+bsxfun(@plus,(0:NTRC-1).*(240+NT*4),uint64([1:4,21:24,37:40,73:76,81:84])');
    H_CSP=reshape(typecast(buff(i(:)),'int32'),5,NTRC);
    % Swap Endianess
    if E;
        H_CSP=swapbytes(H_CSP);
    end;
    H_CSP=double(H_CSP);
    H_CSP(4:5,:)=H_CSP(4:5,:)/1000;
    H_CSP=H_CSP([1,3,4,5,2],:);
end;
%% File Headers
if nargout>2;
    sx=unique(H_CSP(3,:));
    if(numel(sx)<2);
        geo=[double(dt)/1000,double(nt),numel(sx),max(histc(H_CSP(3,:),sx)),0,H_CSP(4,2)-H_CSP(4,1),H_CSP(3,1)];
    else
        geo=[double(dt)/1000,double(nt),numel(sx),max(histc(H_CSP(3,:),sx)),sx(2)-sx(1),H_CSP(4,2)-H_CSP(4,1),H_CSP(3,1)];
    end;
end;
%% IBM -> IEEE
function d=ibm2ieee(ibmf)
% Name:         ibm2ieee
% Abstract:     convert a matrix of IBM/360 32-bit floats
%               to IEEE doubles.
%
%               IBMF is the matrix of IBM/360 32-bit floats each
%               stored as a 32 bit unsigned big-endian integer
%               in a MATLAB double.
%
%               The format of a IBM/360 32-bit float is:
%                  sign 7-bit exponent  24 bit fraction
%                  The base is 16. The decimal point is to
%                  the left of the fraction. The exponent is
%                  biased by +64.
%
%               The basic idea is that you use floating point on
%               the various fields.
%
%               ieee = sign * 16 ^ (exponent - 64) * fraction / 16 ^ 6
%
% By:           Martin Knapp-Cordes
%               The MathWorks, Inc.
%
% Date(s):      Jun 95 - 28, 29

% $Revision: 1.1 $  $Date: 2002/05/28 17:12:55 $
% $Id: ibm2ieee.m,v 1.1 2002/05/28 17:12:55 gary Exp $
% $Log: ibm2ieee.m,v $
% Revision 1.1  2002/05/28 17:12:55  gary
% Initial import into CVS
%
% Revision 3.0  2000/06/13 19:20:32  gilles
% Release 3
%
% Revision 2.0  1999/05/21 18:45:47  mah
% Release 2
%
% Revision 1.1  1999/01/11 19:14:47  kay
% Initial revision
%
%----------------------------------------------------------------------------
%
if (nargin ~= 1)
	error ('Wrong number of arguments.');
elseif (isempty(ibmf))
    error ('Argument is an empty matrix');
end
aibmf=sprintf('%08x',ibmf);
hexd=sscanf(aibmf,'%1x%1x%6x',[3,inf]);
d=(1-(hexd(1,:)>=8).*2).*16.^((hexd(1,:)-(hexd(1,:)>=8).*8).*16+hexd(2,:)-70).*hexd(3,:);
d=reshape((1-(hexd(1,:)>=8).*2).*16.^((hexd(1,:)-(hexd(1,:)>=8).*8).*16+hexd(2,:)-70).*hexd(3,:),size(ibmf));
end
end