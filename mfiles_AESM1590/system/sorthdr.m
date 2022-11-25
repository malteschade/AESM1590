function H_SRT = sorthdr(H_SHT, sortkey1, sortkey2)
% [H_SRT] = SORTHDR(H, sortkey1, sortkey2) 
%
% SORTHDR sorts only  the header of a dataset according to the value of 
% 'sortkey1' as primary sort, and within 'sortkey1' the header is sorted 
% by the value of 'sortkey2'. (sortkey1 must be different from sortkey2).
%
% Valid values for sortkey are:
%    2 = Common Offset
%    3 = Common Shot
%    4 = Common Receiver
%    5 = Common MidPoint (CMP)
%
% Input:        H - header to be sorted
%        sortkey1 - primary sorkey
%        sortkey2 - secondary sortkey (optional)
% output:   H_SRT - sorted header
%
% See also SORTDATA.

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

% transpose, sort and transpose back
Transposed=H_SHT';		% Transpose the headered data matrix,
                        % because sorting can only be done on rows.
SrtdTranspd=sortrows(Transposed,[sortkey1 sortkey2]);
                        % Sort on headerrow. 
H_SRT=SrtdTranspd';		% Transpose the sorted data matrix to its
                        % original position.

