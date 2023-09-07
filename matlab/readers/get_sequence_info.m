function [ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( fasta_file );
% [ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( fasta_file );
% Read in sequences and Eterna information
%
% Inputs
%  fasta_file = Fasta file describing all sequences in library, with tab
%  delimited headers containing numerical id, title, author.
%
% Outputs
%  ids = (list of numbers) Eterna ids (column 1 in header)
%  titles = (cell of strings) titles (column 2 in header)
%  authors = (cell of strings)  Eterna authors (column 3 in header)
%  headers = (cell of strings)  Full headers
%  sequences = (cell of strings) RNA sequences
%  id_strings = (cell of strings) ids as strings
%
% TODO: make this a bit more robust -- should be OK if FASTA file is not
% formatted correctly.
%
% (C) R. Das, HHMI/Stanford University 2023.

fprintf( 'Reading sequence FASTA file...: %s.\n', fasta_file)
seqs = fastaread(fasta_file);

ids = [];
id_strings = repmat({''},1,length(seqs));
titles = repmat({''},1,length(seqs));
authors = repmat({''},1,length(seqs));
for i =1:length(seqs)
    headers{i} = strrep(seqs(i).Header,'%23','#');
    sequences{i} = strrep(seqs(i).Sequence,'T','U');
    cols = strsplit(seqs(i).Header,'\t');
    id_strings{i} = cols{1};
%     ids(i) = NaN;
%     id = str2num(cols{1});
%     if ~isempty(id) ids(i) = id; end;
    if length(cols) > 1; titles{i}=strip(strrep(cols{2},'%23','#')); end
    if length(cols) > 2; authors{i}=strip(strrep(cols{3},'%23','#')); end;
end
ids = str2double(id_strings);
if length(unique(titles))<length(unique(authors)) % swap
    x = titles;
    titles = authors;
    authors = x;
end
fprintf( 'Read in %d sequences from %s.\n', length(sequences),fasta_file)