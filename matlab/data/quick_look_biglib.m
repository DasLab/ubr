function d = quick_look_ubr_biglib(filedir,sequence_file,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,options, tags );
% d = quick_look_ubr_biglib(filedir,sequence_file,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,options, tags );
%
%  Wrapper around quick_look_ubr for big libraries (1M or more sequences).
%  TODO: instead of concatenating in the script, just output chunks to HDF5
%    file
% 
% Inputs
%  filedir       = Directory name with merged UBR output (files like
%       RTB002.muts.txt, RTB002.coverage.txt, etc.)
%  sequence_file = FASTA file describing all sequences in library. Ideally with tab
%       delimited headers containing id, title, author. [Can be FASTA
%       struct read in by fastaread() instead of file.]
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
%  options    = cell of strings specifying possible options. (Default: {})
%                Current option strings are:
%                 'output_all': return all the intermediate data arrays 
%                    in output struct 
%                 'focus_on_shape_nomod': do not use counts across all
%                     conditions to normalize, just the ones specified in
%                     shape_nomod_idx
%                 'no_spread_deletions': for same-nucleotide stretches,
%                   deletions are ambigous and here, by default, are spread across the
%                   stretch according to the mutation signal. Supply
%                   'no_spread_deletions' to suppress this procedure.
%                 'no_raw_counts': do not read in raw_counts or use to
%                    spread deletions
%                 'no_figures': Don't make figures.
%                 'no_print': Don't export figures, which can take
%                       a long time
%                 'no_norm': Don't apply normalization
%                 'no_GA': Exclude G to A mutations
%                 'only_GA': Only count G to A mutations (e.g. DMS at N7)
%                 'only_AG': Only count A to G mutations (e.g. deaminase)
% tags         = cell of strings: tags to use instead of names of hdf5 or
%                      .txt.gz files
%
% Output
%  d = MATLAB struct with the following fields:
%   .sequences = (cell of Ndesigns strings) RNA sequences
%   .r_norm     = [Ndesigns x Nres x Nconditions] Reactivity matrix, normalized.
%   .r_norm_err = [Ndesigns x Nres x Nconditions] Reactivity error matrix, normalized.
%   .conditions = {cell of Nconditions strings} tags for each modificaiton condition in r.
%                   Derived from tags, which have been inferred from files in directory,
%                   and user inputted shape_nomod_idx.
%   .reads    = [Ndesigns x Nconditions] Total reads of signal and background 
%                 contributing to each design under each condition
%   .signal_to_noise = [Ndesigns x Nmodconditions] Estimated signal-to-noise for each
%             design under each experimental modification condition. [Mean of signal
%             divided by mean of error across design]
%   .BLANK_OUT5 = Ignore this number of 5' residues in later analysis and output.
%   .BLANK_OUT3 = Ignore this number of 3' residues in later analysis and output.   
%   .ids = (list of numbers) Eterna ids (column 1 in header)
%   .titles = (cell of Ndesigns strings) titles (column 2 in header)
%   .authors = (cell of Ndesigns strings)  Eterna authors (column 3 in header)
%   .headers = (cell of Ndesigns strings)  Full headers
%   .id_strings = (cell of Ndesigns strings) ids as strings
%   .shape_nomod_idx = {cell of Nconditions pairs of indices} Input shape_nomod_idx
%   .tags = {cell of Ntags strings} could have been inferred from files in directory
%   .coverage = [Ndesigns x Ntags] Total reads for each design under each
%               experimental condition
%   .total_coverage = [Ndesigns] Total reads for each design summed across each
%               experimental condition
%   .structures = (cell of Ndesigns strings) holding structures. Empty if 
%                     no structure_csv_file was specified.
%   .structure_map = [Ndesigns x Nres] matrix of 0,1 for paired/unpaired
%                       across structures. Empty if 
%                     no structure_csv_file was specified.
%   .filedir       = (string) Input filedir
%   .sequence_file = (string) Input sequence_file
%   .structure_csv_file = (string) Input structure_csv_file ('' if not specified)
%
% (C) Rhiju Das, Stanford University & HHMI, 2025.
%
if length(filedir)==0; filedir= './'; end;
if ~exist( 'BLANK_OUT5','var') BLANK_OUT5 = []; end;
if ~exist( 'BLANK_OUT3','var') BLANK_OUT3 = []; end;
if ~exist( 'structure_csv_file','var') structure_csv_file = ''; end;
if ~exist( 'shape_nomod_idx','var') | isempty(shape_nomod_idx); shape_nomod_idx = {}; end;
if ~iscell(shape_nomod_idx) & isnumeric(shape_nomod_idx); shape_nomod_idx = cell(shape_nomod_idx); end;
if ~exist( 'options', 'var') options = {}; end;
if ~exist('tags','var') tags = []; end; 

if isstruct(sequence_file)
    fasta = sequence_file;
else
    if ~exist(sequence_file,'file'); sequence_file = [filedir,'/',sequence_file]; end
    assert(exist(sequence_file,'file'));
    fprintf('Reading in FASTA file %s...\n',sequence_file)
    fasta = fastaread( sequence_file );
end

CHUNK_SEP = 1e4; chunk_size = 1e4; % full run.
Ntot = length(fasta); Nseq = Ntot;
num_chunks = ceil(Ntot/CHUNK_SEP);
for q = 1:num_chunks;
    chunk_start = 1+(q-1)*CHUNK_SEP;
    chunk_end = min(chunk_start + chunk_size-1,Nseq);
    seq_range = [chunk_start,chunk_end];
    all_d{q} = quick_look_ubr(filedir,fasta,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options,seq_range,tags);
    fprintf('\nCOMPLETED: %d of %d\n\n',q,num_chunks)
end

% Concatenate everything into one data structure!
d = concat_ubr_data( all_d );
output_ubr_stats_summary( d );




