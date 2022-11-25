function [sorted, H_SRT] = sortdata(data, H_SHT, sortkey1, sortkey2)
% [data_SRT, H_SRT] = SORTDATA(data, H, sortkey1, sortkey2) 
%
% SORTDATA sorts data and headers H of a input dataset (e.g. read in
% with SEGYREAD) according to the value of 'sortkey1' as primary sort. 
% Within 'sortkey1' the header is sorted by the value of 'sortkey2'. 
% (sortkey1 must be different from sortkey2). 
% 
% Valid values for sortkey are:
%    2 = Common Offset
%    3 = Common Shot
%    4 = Common Receiver
%    5 = Common MidPoint (CMP)
%
% Input:         data - dataset to be sorted
%                   H - its header
%            sortkey1 - primary sorkey
%            sortkey2 - secondary sortkey (optional)
% Output:    data_SRT - sorted dataset
%               H_SRT - its headers
%
% See also PLOTHDR, SEGYREAD

% fill in some defaults if not specified
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

% Display some text
disp(' ')
disp('  SORT is sorting data ... ')
% disp('  This may take some time (depending on the amount of available memory) ...') 

% Read the amount of time-samples and traces from the size of the datamatrix
[nt,ntr]=size(data);
 
% Put the 5 headerrows on top of the data
B=[H_SHT; data];

% Transpose, sort and transpose back
A=B';  	        	% Transpose the headered data matrix,
                    % because sorting can only be done on rows.
B=sortrows(A,[sortkey1 sortkey2]);
                    % Sort on headerrow. 
A=B';   		    % Transpose the sorted data matrix to its
                    % original position.

% Output the sorted headers
H_SRT(1:5,:)=A(1:5,:);

% Strip the headers from the sorted data
sorted(1:nt,:)=A(6:(nt+5),:);

% Display some text
% disp(' ')
disp('  Sorting completed!')
% disp('  Remark: use the CLEAR VARIABLENAME command to clear the')
% disp('  shot-sorted input dataset to save memory (recommended).')
% disp('  Use READSEGY if you need that dataset again.')
disp(' ')
