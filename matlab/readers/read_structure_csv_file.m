function [ structures, structure_map ] = read_structure_csv_file( structure_csv_file, sequences  );
% [ structures, structure_map ] = read_structure_csv_file( structure_csv_file, sequences  );
%
% Inputs
%  structure_csv_file = CSV file holding coumns for 'structure', and
%                         'sequence'
%  sequences = cell of strings of sequences, to check that sequences match
%                    up!
%
% Outputs
%  structures = cell of strings holding structures.
%  structure_map = [Ndesigns x Nres] matrix of 0,1 for paired/unpaired
%                       across structures
%
% (C) R. Das, HHMI/Stanford University 2023.

structure_map = [];
structures = {};
if isempty(structure_csv_file); return; end;
x = readtable(structure_csv_file);
for n=1:height(x)
    i=find(strcmp(sequences,x.sequence{n}));
    assert(length(i)==1);
    if exist('sanitize_structure','file')
        structures{i} = sanitize_structure( x.structure{i} );
    else
        structures{i} = x.structure{i};
    end
    structure_map(i,:) = zeros(length(structures{i}),1);
    structure_map(i,strfind(structures{i},'.')) = 1;
end
fprintf( 'Read in %d structures from %s.\n', length(structures),structure_csv_file);
if ~exist('sanitize_structure','file')
    fprintf('Note that structures were *not* sanitized. If you want sanitize_structure to be used, please install the OpenKnotScore repository and re-run.\n')
end
