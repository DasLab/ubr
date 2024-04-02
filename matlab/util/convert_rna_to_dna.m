function convert_rna_to_dna(fastafile)
% convert_rna_to_dna(fastafile)
%  Converts U to T and prepend T7 promoter ('TTCTAATACGACTCACTATA')
%  Name of output file will have prefix of '.DNA.fa' or '.DNA.fasta'
%
% Input:
%  fastafile = name of fasta file, must end in '.fa' or 'fasta'
%
% (C) R. Das, Stanford & HHMI, 2024

T7promoter = 'TTCTAATACGACTCACTATA';
fprintf('Reading file: %s\n',fastafile);
fasta = fastaread(fastafile);

fasta_new = fasta;
for i = 1:length(fasta)
    fasta_new(i).Sequence = [T7promoter,strrep(fasta(i).Sequence,'U','T')];
end

outfile = strrep(fastafile,'.fa','.DNA.fa');
fprintf('Writing %d lines to %s\n',length(fasta),outfile)
tic
fastawrite(outfile,fasta_new);
toc