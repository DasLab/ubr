function d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options,seq_range);
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
%                 'no_GA': Exclude G to A mutations
%                 'only_GA': Only count G to A mutations
% seq_range    =  [start_idx, end_idx] two integers with range of sequences
%                   to read in. If oneinteger is specified, assumed to be maximum number of 
%                   sequences to read in. Default [], read all sequences.
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
[m,c,rc,tags] = read_ubr_output( filedir,[],use_raw_counts, 0, seq_range);
if isempty(m); finish_quick_look(); return; end;

shape_nomod_idx = update_shape_nomod_idx( shape_nomod_idx, tags);
if isempty(shape_nomod_idx); finish_quick_look(); return; end;

focus_on_shape_nomod = any(strcmp(options,'focus_on_shape_nomod'));
if focus_on_shape_nomod; [m,c,rc,tags,shape_nomod_idx] = focus_on_idx(m,c,rc,tags,shape_nomod_idx); end;

[ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( sequence_file, seq_range);
[ structures, structure_map ] = read_structure_csv_file( structure_csv_file, sequences );
[BLANK_OUT5, BLANK_OUT3] = figure_out_BLANK_OUT( BLANK_OUT5, BLANK_OUT3, sequences );

spread_deletions = ~any(strcmp(options,'no_spread_deletions'));
no_GA = any(strcmp(options,'no_GA'));
only_GA = any(strcmp(options,'only_GA'));
if ~use_raw_counts; rc = m; end;
if spread_deletions 
    % default is to spread deletions.
    [r,r_err,f,f_err,coverage,signal_to_noise] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,sequences,no_GA,only_GA);
else
    [r,r_err,f,f_err,coverage,signal_to_noise] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,[],no_GA,only_GA);
end

for i = 1:length(shape_nomod_idx); conditions{i} = tags{ shape_nomod_idx{i}(1)}; end

total_coverage = sum(coverage,2);
norm_idx = figure_out_idx_for_normalization( total_coverage );
r_norm = []; r_norm_err = [];
if length(norm_idx)>0; 
    [r_norm, r_norm_err,~,norm_val] = normalize_reactivity(r,r_err,norm_idx,BLANK_OUT5, BLANK_OUT3, conditions, [], sequences );
end

for i = 1:length(shape_nomod_idx)
    reads(:,i) = sum(coverage(:,shape_nomod_idx{i}),2);
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

tic
fprintf('\nCreating figures (provide no_figures in options to skip)...\n')
Ntags = size(m,3);
toggle_to_figure(1);
set(gcf,'color','white','Position',[100   100  max(50+Ntags*100,800)   526],'name','fraction reacted, first 500 designs')
Nplot = min(size(f,1),500);
for i = 1:Ntags
    subplot(1,Ntags,i)
    imagesc(sum(f(1:Nplot,:,i),3),[0 0.1]);
    colormap(1-gray(100));
    title(strsplit(tags{i},'_'),'interpreter','none')
    ylim([0.5 Nplot+0.5])
    colorbar('Location','SouthOutside');
end

%% Make plot of signal-to-noise vs. coverage
toggle_to_figure(2);
set(gcf,'color','white','position',[450   991   355   322],'name','S/N vs. reads'); clf
clf
for k = 1:size(reads,2)
    semilogx(reads(:,k),signal_to_noise(:,k),'.'); hold on
    xlabel('Reads (signal+background channels)')
    ylabel('Signal/noise');
end
hold off
legend(conditions,'Interpreter','none');

%% Histogram of signal to noise
toggle_to_figure(3);
set(gcf,'color','white','position',[694   960   400   600],'name','Mean S/N')
clf
for i = 1:size(r,3)
    subplot(size(r,3),1,i);
    s2n = signal_to_noise(:,i); s2n( find(isnan(s2n)) ) = 0.0;
    histogram( s2n )
    xlabel( 'Signal/noise' );
    cols = strsplit(what(filedir).path,'/');
    dirname = strjoin(cols(end-1:end),'/');
    title( [conditions{i},' Mean signal/noise = ',num2str(mean(s2n))],'interp','none');
end

%% Make heat map, up to 500 with high signal to noise 
toggle_to_figure(4);
set(gcf,'color','white','name','first designs with good S/N (up to 500)')
clf
good_idx = find( signal_to_noise(:,end)>=1.0 & reads(:,end) > 100);
Nplot = min(length(good_idx),500);
good_idx = good_idx(1:Nplot);
make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, conditions);

%% Make heat map, up to 10000 with high signal to noise
toggle_to_figure(5);
good_idx = find( signal_to_noise(:,end)>=1.0 & reads(:,end) > 100);
if length(good_idx)>500
    set(gcf,'color','white','name','first designs with good S/N (up to 10000)')
    clf
    Nplot = min(length(good_idx),10000);
    good_idx = good_idx(1:Nplot);
    make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, conditions);
else
   close(5);
end


%% Take a close look at one of the constructs with high apparent signal to noise
toggle_to_figure(6);
set(gcf,'color','white','position',[599   477   560   420],'name','Top S/N design')
clf
[~,idx] = max(sum(signal_to_noise,2));
for i = 1:length(shape_nomod_idx)
    subplot(length(shape_nomod_idx),1,i);
    cidx = shape_nomod_idx{i}; %(end:-1:1);
    plot( squeeze(f(idx,:,cidx)),'LineWidth',2 )
    h = legend( tags{cidx});
    set(h,'interp','none');
end
h=title(headers(idx));
set(h,'interpreter','none')

%% Look through each of the conditions - mutational profiles (mean over designs)
if use_raw_counts
    show_mut_type_analysis(mut_rate_matrix,rfcount_mut_rate_profiles,tags,tags,BLANK_OUT5, BLANK_OUT3);
end
toc

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



