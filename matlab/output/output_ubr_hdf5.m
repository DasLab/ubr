function output_ubr_hdf5( filename, d, VERBOSE );
% output_ubr_hdf5( filename, d );
%
%  Works generally on any MATLAB workspace, but
%   can have special treatment of shape_nomod_idx and
%   other variables in  quick_look_ubr().
%
% Inputs:
%  filename = name of output HDF5-formatted file.
%  d        = MATLAB workspace 
%  VERBOSE  = verbose output (Default 0)
%
% (C) R. Das, Stanford University & HHMI 2024
if ~exist('VERBOSE','var'); VERBOSE = 0; end; 
fields = fieldnames(d);
if exist(filename,'file'); delete(filename); end;
if VERBOSE; fprintf('Outputting to HDF5 file: %s\n',filename); end;
for i = 1:length(fields)
    fieldName = fields{i};
    data = d.(fieldName);
    if strcmp(fieldName,'shape_nomod_idx'); data = cell2mat(data); end;
    if ~isempty(data)
        if VERBOSE; fprintf('Outputting %s\n',fieldName); end;
        if iscell(data); data = string(data); end;
        if isnumeric(data) | isstring(data)
            sz = fliplr(size(data));
            h5create(filename,['/',fieldName],sz,'Datatype',class(data));
            h5write(filename,['/',fieldName],permute(data,length(sz):-1:1) );
        elseif ischar(data)
            data = string(data);
            sz = size(data);
            h5create(filename,['/',fieldName],sz,'Datatype','string');
            h5write(filename,['/',fieldName], data );
        else
            warning( ['Did not output non-empty field: ',fieldName] );
        end
    else
        %warning( ['Did not output empty field: ',fieldName] );
    end
end
