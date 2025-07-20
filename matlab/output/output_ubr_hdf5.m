function output_ubr_hdf5( filename, d, VERBOSE, fields, section );
% output_ubr_hdf5( filename, d );
%
%  Works generally on any MATLAB struct, but
%   can have special treatment of shape_nomod_idx and
%   other variables in  quick_look_ubr().
%
% Inputs:
%  filename = name of output HDF5-formatted file.
%  d        = MATLAB struct (or optionally, cell of struts)
%
% Optional:
%  VERBOSE  = verbose output (Default 1)
%  fields   = cell of strings: which fields to use (Default: all fields).
%               Can be cell of cell of strings if multiple structs are
%               input for d
%  section = (string) names of section(s) (Default: [], no sections). 
%                Can be cell of strings if multiple structs are input for d.
%
% (C) R. Das, Stanford University & HHMI 2024-2025
if ~exist('VERBOSE','var'); VERBOSE = 1; end; 
if ~exist('fields','var'); fields = []; end; % all fields
if ~exist('section','var'); section = []; end; % no section

all_d = d; all_fields = fields; all_sections = section;
if isstruct(d); all_d = {all_d}; all_fields = {all_fields}; all_sections = {section}; end;

if exist(filename,'file'); delete(filename); end;
if VERBOSE; fprintf('Outputting to HDF5 file: %s\n',filename); end;

for i = 1:length(all_d);
    d = all_d{i};
    fields = all_fields{i};
    section = all_sections{i};
    output_struct( d, fields, section, filename, VERBOSE )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_struct( d, fields, section, filename, VERBOSE )
if isempty(fields)
    fields = fieldnames(d);
end
section_tag = '';
if ~isempty(section); 
    section_tag = ['/',section]; 
    fprintf('\nOutputting section: %s\n',section)
end
for i = 1:length(fields)
    fieldName = fields{i};
    data = d.(fieldName);
    if strcmp(fieldName,'shape_nomod_idx'); data = cell2mat(data); end;
    if ~isempty(data)
        if VERBOSE; fprintf('Outputting %s\n',fieldName); end;
        if iscell(data); data = string(data); end;
        if isnumeric(data) | isstring(data)
            sz = fliplr(size(data));
            h5create(filename,[section_tag,'/',fieldName],sz,'Datatype',class(data));
            h5write(filename,[section_tag,'/',fieldName],permute(data,length(sz):-1:1) );
        elseif ischar(data)
            data = string(data);
            sz = size(data);
            h5create(filename,[section_tag,'/',fieldName],sz,'Datatype','string');
            h5write(filename,[section_tag,'/',fieldName], data );
        else
            warning( ['Did not output non-empty field: ',fieldName] );
        end
    else
        %warning( ['Did not output empty field: ',fieldName] );
    end
end
