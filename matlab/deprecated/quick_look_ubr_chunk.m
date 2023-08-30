function d_chunks = quick_look_ubr_chunk(filedir,sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options,chunk_size);
% Read in data from ubr (ultraplex-bowtie2-RNAframework pipeline),
%  pre-chunk sequences and data files.
%
% Usage:
%  d_chunks = quick_look_ubr_chunk(filedir,sequence_file,shape_nomod_idx);
%  d_chunks = quick_look_ubr_chunk(filedir,sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options,chunk_size);
%
%
% Inputs
%  filedir       = Directory name with merged UBR output (files like
%       RTB002.muts.txt, RTB002.coverage.txt, etc.)
%  sequence_file = FASTA file describing all sequences in library. Ideally with tab
%       delimited headers containing id, title, author.
%  shape_nomod_idx = cell of Nconditions pairs of tags or indices into experimental conditions
%            in order [shape_i,nomod_i]. 
%                Example (tags): {'RTB028_Agilent_15k_DMS','RTB029_Agilent_15k_MaraNomod'},{'RTB030_Agilent_15k_2A3','RTB031_Agilent_15k_SSIINomod'}
%                Example (indices): {[1,2],[3,4]}.
%            These numbers refer to the order in which the tags for the
%            files appear in filedir. Within each pair of tags/indices, put
%            the modification condition first (e.g., SHAPE or DMS), and the
%            no modification condition second.
%
% Optional inputs
%  structure_csv_file = structure_csv_file = CSV file holding columns for 
%         'structure' and 'sequence'. Provide as '' or [] if not available.
%  BLANK_OUT5 = Ignore this number of 5' residues in signal_to_noise
%                  calculation. Default [], which triggers to 26 if using
%                  sequences that all begin with GGGAACGACUCGAGUAGAGUCGAAAA, 0 otherwise.                
%  BLANK_OUT3 = Ignore this number of 3' residues in signal_to_noise
%                  calculation.  Default [], which triggers to 26 if using
%                  sequences that all end with AAAGAAACAACAACAACAAC, 0 otherwise.    
%  options     = cell of strings specifying possible options. (Default: {})
%                Current option strings are:
%                 'output_all': return all the intermediate data arrays 
%                    in output struct 
%                 'no_spread_deletions': for same-nucleotide stretches,
%                   deletions are ambigous and here, by default, are spread across the
%                   stretch according to the mutation signal. Supply
%                   'no_spread_deletions' to suppress this procedure.
%                 'no_raw_counts': do not read in raw_counts or use to
%                    spread deletions
%                 'no_figures': Don't make figures.
%                 'no_print': Don't export figures, which can take
%                       a long time
%  chunk_size    = [integer] number of sequences per chunk (Default 10000)
%
% Output
%  d_chunks = MATLAB struct with the following fields:
%   .all_chunk_idx = cell of num_chunks arrays, each of which contains
%                    integer indices (chunk_idx) of which sequences ended up
%                    in which chunk.
%   .workspace_files = .MAT formatted save files for each individual chunk
%                    containing single struct "d" with r_norm, sequences,
%                    etc. for each individual chunk.
%   .sequences = (cell of Ndesigns strings) RNA sequences
%   .ids = (list of numbers) Eterna ids (column 1 in header)
%   .titles = (cell of Ndesigns strings) titles (column 2 in header)
%   .authors = (cell of Ndesigns strings)  Eterna authors (column 3 in header)
%   .headers = (cell of Ndesigns strings)  Full headers
%   .id_strings = (cell of Ndesigns strings) ids as strings
%
%  If "output_all" specified in options, d will also contain these potentially big
%        arrays:
%   .m  = [Ndesigns x Nres x Ntags] mutation counts 
%   .c  = [Ndesigns x Nres x Ntags] coverage counts 
%   .rc = [Ndesigns x Nres x Nmuttypes x Ntags] coverage counts, divided by
%            mutation type (comes out of RNA-framework with -orc flag)
%   .r     = [Ndesigns x Nres x Nmodconditions] Reactivity matrix, as fraction at each position that leads to mutation, background subtracted.
%   .r_err = [Ndesigns x Nres x Nmodconditions] Error of reactivity matrix, as fraction at each position that leads to mutation.
%   .f     = [Ndesigns x Nres x Ntags] Fraction at each position that leads to mutation.
%   .f_err = [Ndesigns x Nres x Ntags] Error ofraction at each position that leads to mutation.
%
% (C) R. Das, HHMI/Stanford University 2023.
if ~exist( 'BLANK_OUT5','var') BLANK_OUT5 = []; end;
if ~exist( 'BLANK_OUT3','var') BLANK_OUT3 = []; end;
if ~exist( 'structure_csv_file','var') structure_csv_file = ''; end;
if ~exist( 'options', 'var') options = {}; end;
if ~exist( 'chunk_size','var'); chunk_size = 10000; end;

d_chunks = struct();
assert(exist(filedir,'dir'));
assert(exist(sequence_file,'file'));
if ~isempty(structure_csv_file) assert(exist(structure_csv_file,'file')); end;

all_chunk_idx = ubr_chunk( filedir, sequence_file, chunk_size );
num_chunks = length(all_chunk_idx);

%% Read in the data and save...
options = [options,{'no_print'}];
for i = 1:num_chunks
    chunk_dir = sprintf('%s/%04d/','CHUNKS',i);
    [~,basename,ext]=fileparts(sequence_file);
    chunk_sequence_file = sprintf('%s/%s%s',chunk_dir,basename,ext);
    d = quick_look_ubr(chunk_dir,chunk_sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options);
    workspace_files{i} = sprintf('workspace_%04d.mat',i);
    fprintf('\nSaving to %s...\n',workspace_files{i})
    save( workspace_files{i},'d' );
end

%% Read in sequence information
[ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( sequence_file );

d_chunks.all_chunk_idx = all_chunk_idx;
d_chunks.sequences = sequences;
d_chunks.workspace_files = workspace_files;
d_chunks.ids = ids;
d_chunks.titles = titles;
d_chunks.authors = authors;
d_chunks.headers = headers;
d_chunks.id_strings = id_strings;



