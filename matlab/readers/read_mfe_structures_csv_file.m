function [x,structure_tags, structure_sets,structure_map] = read_mfe_structures_csv_file( structure_sets_csv_file, ordered_sequences )
% [x,structure_tags, structure_sets,structure_map] = read_structure_sets_csv_file( structure_sets_csv_file, ordered_sequences );
%
% MOVED TO OpenKnotScore!
%
% Inputs
%  structure_sets_csv_file = csv file with columns like "_mfe" holding
%                        structure predictions from different packages in dot bracket notation.
%  ordered_sequences = [Optional] list of sequences. If provided, structures read
%          in from .csv file will be reordered based on sequence column to match ordering in
%          ordered_sequence
%
% Outputs
% x = csv file in MATLAB Table object
% structure_tags = tags for each structure/mfe column
% structure_sets = cell of cell of strings of predicted structures
% structure_map = [Ndesign x Nres x Npackage] matrix of 0,1 for
%               paired/unpaired in each package structure prediction
%
% (C) R. Das, HHMI/Stanford University 2023.

if ~exist('ordered_sequences','var');
  ordered_sequences = [];
end;

if ~exist( 'read_structure_sets_csv_file', 'file' )
  fprintf( 'ERROR! You need to install OpenKnotScore from https://github.com/eternagame/OpenKnotScore, and add OpenKnotScore/matlab and subdirectories to your path.' );
end

[x,structure_tags, structure_sets,structure_map] = read_structure_sets_csv_file( structure_sets_csv_file, ordered_sequences );

fprintf( 'WARNING! read_mfe_structures_csv_file has been deprecated in favor of read_structure_sets_csv_file in the OpenKnotScore repository at https://github.com/eternagame/OpenKnotScore\n' );

