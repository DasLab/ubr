function d = read_ubr_hdf5( filename, fields, seq_range );
% d = read_ubr_hdf5( filename, fields, seq_range );
%
% Reads HDF5 file into MATLAB workspace.
%
% Inputs:
%  filename = name of  HDF5-formatted file.
%  fields = [optional] [cell of strings] names of variables to read in.
%               (Default: all fields)
%  seq_range = [optional] [start_idx,end_idx] to read in . (Default: all)
%
% Output:
%  d        = MATLAB workspace 
%
% (C) R. Das, Stanford University & HHMI 2024
if ~exist('fields','var') fields = []; end;
if ischar(fields); fields = {fields}; end;
if ~exist( 'seq_range', 'var') seq_range = []; end;
if length(seq_range)==1; seq_range = [1 seq_range(1)]; end;

% main loops.
d = read_ubr_hdf5_data( filename, '', fields, seq_range); 

groups = h5info(filename).Groups;
for i = 1:length(groups)
    name = groups(i).Name(2:end); % get rid of '/' at beginning   
    d.(name) = read_ubr_hdf5_data( filename,name, fields, seq_range);
end






