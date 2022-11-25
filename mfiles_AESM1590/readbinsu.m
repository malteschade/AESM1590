function  readbin = readbinsu(fname,ntr,nt)
% readbinsu: this function reads a binary file that has been created from 
% a Seismic-Unix (su) file by stripping the headers
%    fname - name of the bin file
%    ntr   - number of traces in the bin file
%    nt    - number of time samples in the bin file

% Reading the binary file into Matlab.
fid  = fopen(fname,'r');
temp = fread(fid,ntr*nt,'float32',0,'ieee-le');
fclose(fid);
readbin = reshape(temp,nt,ntr);
end