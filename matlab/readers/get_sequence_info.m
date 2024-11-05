function [ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( fasta_file, seq_range );
% [ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( fasta_file, seq_range );
% Read in sequences and Eterna information
%
% Inputs
%  fasta_file = Fasta file describing all sequences in library, with tab
%                 delimited headers containing numerical id, title, author.
%                 Can also be FASTA struct() read in by fastaread().
%  MAX_SIZE   = (integer) Maximum number of sequences to read in [default 0, no max]
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

if ~exist( 'seq_range', 'var') seq_range = []; end;
if length(seq_range)==1; seq_range = [1 seq_range(1)]; end;

if ischar(fasta_file)
    fprintf( 'Reading sequence FASTA file...: %s.\n', fasta_file)
    seqs = fastaread(fasta_file);
else
    assert(isstruct(fasta_file));
    seqs = fasta_file;
end

if ~isempty(seq_range) > 0
    if length(seqs) < seq_range(2); seq_range(2) = length(seqs); end;
    seqs = seqs( [seq_range(1):seq_range(2)] );
end

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
if ischar(fasta_file); fprintf( 'Read in %d sequences from %s.\n', length(sequences),fasta_file); end;

