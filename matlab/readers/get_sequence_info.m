function [ids,titles,authors,headers,sequences] = get_sequence_info( fasta_file );
% [ids,titles,authors,headers,sequences] = get_sequence_info( fasta_file );
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
%
% TODO: make this a bit more robust -- should be OK if FASTA file is not
% formatted correctly.
%
% (C) R. Das, HHMI/Stanford University 2023.

seqs = fastaread(fasta_file);
ids = [];
for i =1:length(seqs)
    headers{i} = strrep(seqs(i).Header,'%23','#');
    sequences{i} = seqs(i).Sequence;
    cols = strsplit(seqs(i).Header,'\t');
    ids(i)=str2num(cols{1});
    titles{i}=strip(strrep(cols{2},'%23','#'));
    authors{i}=strip(cols{3});
end

fprintf( 'Read in %d sequences from %s.\n', length(sequences),fasta_file)