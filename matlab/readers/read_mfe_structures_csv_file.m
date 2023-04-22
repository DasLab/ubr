function [x,mfe_tags, mfe_structures,mfe_structure_map] = read_mfe_structures_csv_file( mfe_structures_csv_file );
% [x,mfe_tags, mfe_structures,mfe_structure_map] = read_mfe_structures_csv_file( mfe_structures_csv_file );
%
% Inputs
%  mfe_structures_csv_file = csv file with columns like "_mfe" holding
%                        structure predictions from different packages in dot bracket notation.
% Outputs
% x = csv file in MATLAB Table object
% mfe_tags = tags for each structure/mfe column
% mfe_structures = cell of cell of strings of predicted structures
% mfe_structure_map = [Ndesign x Nres x Npackage] matrix of 0,1 for
%               paired/unpaired in each package structure prediction
%
% (C) R. Das, HHMI/Stanford University 2023.

x = readtable(mfe_structures_csv_file);
mfe_tags = {}; count = 0;
for n = 1:length(x.Properties.VariableNames);
    tag = x.Properties.VariableNames{n};
    structures = table2cell(x(:,n));
    if contains(tag,'_mfe') | ...
            contains(tag,'structure') | ...
            ( contains(tag,'hfold') & ...
            ~contains(tag,'time')  | ...
            strcmp(tag,'eterna_nupack') ) | ...
            (ischar(structures{1}) & contains(structures{1},'('))
        if ~ischar( structures{1}); continue; end;
        count = count + 1;
        mfe_tags{count} = strrep(strrep(tag,'__mfe',''),'_mfe','');
        mfe_structures{count} = structures;
        fprintf( 'Sanitizing %d structures for %s...\n',size(x,1),mfe_tags{count});
        for i = 1:size(x,1) % loop over designs
            mfe_structures{count}{i} = sanitize_structure( mfe_structures{count}{i} );
        end
    end
end

mfe_structure_map = get_mfe_structure_map( mfe_structures );
