function d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options,seq_range,tags);
% Read in data from ubr (ultraplex-bowtie2-RNAframework pipeline)
%
% Usage:
%  d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx);
%  d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options,seq_range);
%
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
%                 'output_rsep_norm': return rsep_norm with mutA/mutC/mutG/mutT/del/ins as
%                    6 separate tracks
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
% seq_range    =  [start_idx, end_idx] two integers with range of sequences
%                   to read in. If oneinteger is specified, assumed to be maximum number of 
%                   sequences to read in. Default [], read all sequences.
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
%   .rsep_norm     = [Ndesigns x Nres x 6 x Nconditions]
%                        Mut(A/C/G/U),del,ins, normalized. Sum of 5 should
%                        equal r_norm.
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
if length(filedir)==0; filedir= './'; end;
if ~exist( 'BLANK_OUT5','var') BLANK_OUT5 = []; end;
if ~exist( 'BLANK_OUT3','var') BLANK_OUT3 = []; end;
if ~exist( 'structure_csv_file','var') structure_csv_file = ''; end;
if ~exist( 'shape_nomod_idx','var') | isempty(shape_nomod_idx); shape_nomod_idx = {}; end;
if ~iscell(shape_nomod_idx) & isnumeric(shape_nomod_idx); shape_nomod_idx = cell(shape_nomod_idx); end;
if ~exist( 'options', 'var') options = {}; end;
if ~exist( 'seq_range', 'var') seq_range = []; end;
if length(seq_range)==1; seq_range = [1 seq_range(1)]; end;

d = struct();
assert(exist(filedir,'dir'));
assert(iscell(shape_nomod_idx));
if ~isstruct(sequence_file)
    if ~exist(sequence_file,'file'); sequence_file = [filedir,'/',sequence_file]; end
    assert(exist(sequence_file,'file'));
end
if ~isempty(structure_csv_file) assert(exist(structure_csv_file,'file')); end;

%% Get logging setup
logfile = sprintf('%s/quick_look_ubr.log',filedir);
if exist(logfile,'file'); delete(logfile); end;
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nPutting output to %s\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',logfile);
diary(logfile); diary on;

%% Data readin
tic
fprintf('Reading and processing data...\n')

use_raw_counts = ~any(strcmp(options,'no_raw_counts'));
[m,c,rc,tags_on_disk] = read_ubr_output( filedir,[],use_raw_counts, 0, seq_range);
if exist('tags','var') & ~isempty(tags); assert( length(tags_on_disk)==length(tags)); else, tags = tags_on_disk; end;

if isempty(m); finish_quick_look(); return; end;

shape_nomod_idx = update_shape_nomod_idx( shape_nomod_idx, tags);
if isempty(shape_nomod_idx); finish_quick_look(); return; end;

focus_on_shape_nomod = any(strcmp(options,'focus_on_shape_nomod'));
if focus_on_shape_nomod; [m,c,rc,tags,shape_nomod_idx] = focus_on_idx(m,c,rc,tags,shape_nomod_idx); end;

[ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( sequence_file, seq_range);
[ structures, structure_map ] = read_structure_csv_file( structure_csv_file, sequences );
[BLANK_OUT5, BLANK_OUT3] = figure_out_BLANK_OUT( BLANK_OUT5, BLANK_OUT3, sequences );

if ~use_raw_counts; rc = m; end;
[r,r_err,f,f_err,coverage,signal_to_noise,~,rsub] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,sequences,options);

for i = 1:length(shape_nomod_idx); conditions{i} = tags{ shape_nomod_idx{i}(1)}; end

total_coverage = sum(coverage,2);
r_norm = r; r_norm_err = r_err; norm_val = ones(1,length(conditions)); norm_idx = [];
if ~any(strcmp(options,'no_norm')); 
    norm_idx = figure_out_idx_for_normalization( total_coverage );
    if length(norm_idx)>0
        [r_norm, r_norm_err,~,norm_val] = normalize_reactivity(r,r_err,norm_idx,BLANK_OUT5, BLANK_OUT3, conditions, [], sequences );
    end
end

for i = 1:length(shape_nomod_idx)
    reads(:,i) = sum(coverage(:,shape_nomod_idx{i}),2);
end

if any(strcmp(options,'output_rsep_norm'))
    mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};
    rsep(:,:,1,:) = sum(rsub(:,:,[4 7 10],:),3); % mut to A
    rsep(:,:,2,:) = sum(rsub(:,:,[1 8 11],:),3); % mut to C
    rsep(:,:,3,:) = sum(rsub(:,:,[2 5 12],:),3); % mut to G
    rsep(:,:,4,:) = sum(rsub(:,:,[3 6  9],:),3); % mut to T
    rsep(:,:,5,:) = sum(rsub(:,:,[14],:),3); % del
    rsep(:,:,6,:) = sum(rsub(:,:,[13],:),3); % ins
    for i = 1:length(norm_val); rsep_norm(:,:,:,i) = rsep(:,:,:,i)/norm_val(i); end;
    rsep_norm(:,[1                 :BLANK_OUT5],:,:) = NaN;
    rsep_norm(:,[(end-BLANK_OUT3+1):end       ],:,:) = NaN;
end

[mut_rate_matrix, rfcount_mut_rate_profiles, coverage_matrix] = get_mut_rate_matrix( m,c,rc );

%% Fill up output struct
d.sequences = sequences;
d.r_norm = r_norm;
d.r_norm_err = r_norm_err;
d.conditions = conditions;
d.reads = reads;
d.signal_to_noise = signal_to_noise;
d.BLANK_OUT5 = BLANK_OUT5;
d.BLANK_OUT3 = BLANK_OUT3;
d.ids = ids;
d.titles = titles;
d.authors = authors;
d.headers = headers;
d.id_strings = id_strings;
d.shape_nomod_idx = shape_nomod_idx;
d.tags = tags;
d.coverage = coverage;
d.total_coverage = total_coverage;
d.norm_idx = norm_idx;
d.norm_val = norm_val;
d.structures = structures;
d.structure_map = structure_map;
d.filedir = filedir;
if ischar(sequence_file); d.sequence_file = sequence_file; end;
d.structure_csv_file = structure_csv_file;
d.mut_rate_matrix = mut_rate_matrix;
d.rfcount_mut_rate_profiles = rfcount_mut_rate_profiles;
d.coverage_matrix = coverage_matrix;
if ~isempty(find(strcmp(options,'output_rsep_norm'))); d.rsep_norm = rsep_norm; end
if ~isempty(find(strcmp(options,'output_all')));
    d.m = m;
    d.c = c;
    d.rc = rc;
    d.r = r;
    d.r_err = r_err;
    d.f = f;
    d.f_err = f_err;
end

%% Stats summary
output_ubr_stats_summary(d);

%% Make some heatmaps
if any(strcmp(options,'no_figures')); finish_quick_look(); return; end;

make_quick_look_figures( d, f );

% Print out .png's.
no_print = any(strcmp(options,'no_print'));
if no_print; finish_quick_look(); return; end;
fprintf('\nPrinting figures (provide no_print in options to skip)...\n')
tic
print_quick_look_figures( filedir );
toc
finish_quick_look()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BLANK_OUT5, BLANK_OUT3] = figure_out_BLANK_OUT( BLANK_OUT5, BLANK_OUT3, sequences );
test_sequence = strrep(sequences{1},'T','U');
if isempty(BLANK_OUT5)
    if length(test_sequence)>26 & strcmp(test_sequence(1:26),'GGGAACGACUCGAGUAGAGUCGAAAA')
        BLANK_OUT5 = 26;
    elseif length(test_sequence)>46 & strcmp(test_sequence(1:46),'UUCUAAUACGACUCACUAUAGGGAACGACUCGAGUAGAGUCGAAAA')
        BLANK_OUT5 = 46;
    else
        BLANK_OUT5 = 0;
    end
    fprintf( 'Setting BLANK_OUT5 to be %d nucleotides; if this does not look right, re-run with explicit specification of BLANK_OUT5\n',BLANK_OUT5);
end
if isempty(BLANK_OUT3)
    if length(test_sequence)>20 & strcmp(test_sequence(end-19:end),'AAAGAAACAACAACAACAAC')
        BLANK_OUT3 = 20;
    else
        BLANK_OUT3 = 0;
    end
    fprintf( 'Setting BLANK_OUT3 to be %d nucleotides; if this does not look right, re-run with explicit specification of BLANK_OUT3\n',BLANK_OUT3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shape_nomod_idx = update_shape_nomod_idx( shape_nomod_idx, tags);
if isempty(shape_nomod_idx) shape_nomod_idx = num2cell(1:length(tags)); end; % look through all tags.
if length(shape_nomod_idx)==0; return; end;
if iscell( shape_nomod_idx{1} );
    for i = 1:length(shape_nomod_idx);
        if length( shape_nomod_idx{i}) == 0; shape_nomod_idx = [];; end;
        if length( shape_nomod_idx{i}) > 2; shape_nomod_idx = [];; end;
        if ~ischar( shape_nomod_idx{i}{1}); shape_nomod_idx = []; return; end;
        [shape_nomod_idx_new1,ok] = get_idx(shape_nomod_idx{i}{1},tags);
        if ~ok; shape_nomod_idx = []; return; end;
        if length( shape_nomod_idx{i}) == 1; shape_nomod_idx{i} = [shape_nomod_idx_new1]; continue; end;
        [shape_nomod_idx_new2,ok] = get_idx(shape_nomod_idx{i}{2},tags);        
        shape_nomod_idx{i} = [shape_nomod_idx_new1,shape_nomod_idx_new2];
    end
end
for i = 1:length(shape_nomod_idx);
    if length( shape_nomod_idx{i} ) > 1
        fprintf('Condition %d will be %s minus %s.\n', i, tags{shape_nomod_idx{i}(1)}, tags{shape_nomod_idx{i}(2)});
    else
        assert( length(shape_nomod_idx{i}) == 1 );
        fprintf('Condition %d will be %s (no background subtraction).\n', i, tags{shape_nomod_idx{i}} );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx,ok] = get_idx(tag,tags);
idx = find(strcmp(tags,tag));
ok = 1;
if length(idx) ~= 1;
    idx = 0;
    ok = 0;
    fprintf('Could not find shape_nomod_idx tag %s in tags read from disk!\n',tag);
    for i = 1:length(tags)
        fprintf('  %s\n',tags{i});
    end
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,c,rc,tags,shape_nomod_idx] = focus_on_idx(m,c,rc,tags,shape_nomod_idx)

% where are we focusing our attention?
focus_idx = [];
for i = 1:length(shape_nomod_idx)
    focus_idx = union( focus_idx, shape_nomod_idx{i});
end
[focus_idx_sorted, sortidx] = sort(focus_idx);
idx_to_sortidx( focus_idx ) = sortidx;
m = m(:,:,focus_idx);
c = c(:,:,focus_idx);
rc = rc(:,:,:,focus_idx);
tags = tags(focus_idx);
for i = 1:length(shape_nomod_idx)
    shape_nomod_idx{i} = idx_to_sortidx(shape_nomod_idx{i});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function finish_quick_look();
diary off
fprintf('%s\n',datetime);



