function d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options);
% Read in data from ubr (ultraplex-bowtie2-RNAframework pipeline)
%
% Usage:
%  d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx);
%  d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx,structure_csv_file,BLANK_OUT5,BLANK_OUT3,options);
%
%
% Inputs
%  filedir       = Directory name with merged UBR output (files like
%       RTB002.muts.txt, RTB002.coverage.txt, etc.)
%  sequence_file = FASTA file describing all sequences in library. Ideally with tab
%       delimited headers containing id, title, author.
%  shape_nomod_idx = cell of Nconditions pairs of indices into experimental conditions
%            in order [shape_i,nomod_i]. Example: {[3,1],[4,2]}.
%            These numbers refer to the order in which the tags for the
%            files appear in filedir, e.g., if filedir has
%            'RTB002.muts.txt','RTB004.muts.txt','RTB006.muts.txt','RTB008.muts.txt'
%            {[3,1],[4,2]} means to take RTB006 and subtract background from RTB002 for the
%            first condition and take RTB008 and subrtact background RTB004
%            as the second condition. Despite the name 'shape_nomod_idx',
%            the experimental condition can be DMS or other mutational
%            profiling approach besides SHAPE reactions.
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
%                 'no_figures': Don't make figures.
%                 'no_print': Don't export figures, which can take
%                       a long time
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
if ~exist( 'BLANK_OUT5','var') BLANK_OUT5 = []; end;
if ~exist( 'BLANK_OUT3','var') BLANK_OUT3 = []; end;
if ~exist( 'structure_csv_file','var') structure_csv_file = ''; end;
if ~exist( 'options', 'var') options = {}; end;
d = struct();

%% Data readin
tic
fprintf('Reading and processing data...\n')
[m,c,rc,tags] = read_ubr_output( filedir,[],1 );
[ids,titles,authors,headers,sequences,id_strings] = get_sequence_info( sequence_file );
[ structures, structure_map ] = read_structure_csv_file( structure_csv_file, sequences );
[BLANK_OUT5, BLANK_OUT3] = figure_out_BLANK_OUT( BLANK_OUT5, BLANK_OUT3, sequences );
if any(strcmp(options,'no_spread_deletions'))
    % default is to spread deletions.
    [r,r_err,f,f_err,coverage,signal_to_noise] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,sequences);
else
    [r,r_err,f,f_err,coverage,signal_to_noise] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3);
end

for i = 1:length(shape_nomod_idx); conditions{i} = tags{ shape_nomod_idx{i}(1)}; end

total_coverage = sum(coverage,2);
good_idx = figure_out_idx_for_normalization( total_coverage );
r_norm = []; r_norm_err = [];
if length(good_idx)>0; 
    [r_norm, r_norm_err] = normalize_reactivity(r,r_err,good_idx,BLANK_OUT5, BLANK_OUT3, conditions );
end

for i = 1:length(shape_nomod_idx)
    reads(:,i) = sum(coverage(:,shape_nomod_idx{i}),2);
end

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
d.structures = structures;
d.structure_map = structure_map;
d.filedir = filedir;
d.sequence_file = sequence_file;
d.structure_csv_file = structure_csv_file;
if ~isempty(strcmp(options,'output_all'));
    d.m = m;
    d.c = c;
    d.rc = rc;
    d.r = r;
    d.r_err = r_err;
    d.f = f;
    d.f_err = f_err;
end

%% Make some heatmaps
if any(strcmp(options,'no_figures')); return; end;
no_print = ~any(strcmp(options,'no_print'));

tic
fprintf('\n\nCreating figures...\n')
set(figure(1),'color','white','Position',[100   785   554   526],'name','fraction reacted, first 500 designs')
Ntags = size(m,3);
Nplot = min(size(f,1),500);
for i = 1:Ntags
    subplot(1,Ntags,i)
    imagesc(sum(f(1:Nplot,:,i),3),[0 0.1]);
    colormap(1-gray(100));
    title(strsplit(tags{i},'_'),'interpreter','none')
    ylim([0 Nplot])
    colorbar('Location','SouthOutside');
end

%% Make plot of signal-to-noise vs. coverage
set(figure(2),'color','white','position',[450   991   355   322],'name','S/N vs. reads'); clf
clf
for k = 1:size(reads,2)
    semilogx(reads(:,k),signal_to_noise(:,k),'.'); hold on
    xlabel('Reads (signal+background channels)')
    ylabel('Signal/noise');
end
hold off
legend(conditions{k},'Interpreter','none');

%% Histogram of signal to noise
set(figure(3),'color','white','position',[694   960   400   344],'name','Mean S/N')
clf
for i = 1:size(r,3)
    subplot(size(r,3),1,i);
    s2n = signal_to_noise(:,i); s2n( find(isnan(s2n)) ) = 0.0;
    histogram( s2n )
    xlabel( 'Signal/noise' );
    title( {filedir,conditions{i},['Mean signal/noise = ',num2str(mean(s2n))]},'interp','none');
    fprintf( 'Mean signal-to-noise (%s) %f\n',conditions{i},mean(s2n) );
end

%% Make heat map, up to 500 with high signal to noise 
set(figure(4),'color','white','name','first designs with good S/N (up to 500)')
clf
good_idx = find( signal_to_noise(:,end)>=1.0 & reads(:,end) > 100);
Nplot = min(length(good_idx),500);
good_idx = good_idx(1:Nplot);
make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, conditions);

%% Make heat map, up to 10000 with high signal to noise
good_idx = find( signal_to_noise(:,end)>=1.0 & reads(:,end) > 100);
if length(good_idx)>10000
    set(figure(5),'color','white','name','first designs with good S/N (up to 10000)')
    clf
    Nplot = min(length(good_idx),10000);
    good_idx = good_idx(1:Nplot);
    make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, conditions);
end

%% Take a close look at one of the constructs with high apparent signal to noise
set(figure(6),'color','white','position',[599   477   560   420],'name','Top S/N design')
clf
[~,idx] = max(sum(signal_to_noise,2));
for i = 1:length(shape_nomod_idx)
    subplot(length(shape_nomod_idx),1,i);
    cidx = shape_nomod_idx{i}(end:-1:1);
    plot( squeeze(f(idx,:,cidx)),'LineWidth',2 )
    h = legend( conditions{i});
    set(h,'interp','none');
end
h=title(headers(idx));
set(h,'interpreter','none')
toc

%% Look through each of the conditions - mutational profiles (mean over designs)
run_mut_type_analysis( m,c,rc,tags,tags,BLANK_OUT5, BLANK_OUT3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BLANK_OUT5, BLANK_OUT3] = figure_out_BLANK_OUT( BLANK_OUT5, BLANK_OUT3, sequences );
test_sequence = strrep(sequences{1},'T','U');
if isempty(BLANK_OUT5)
    if length(test_sequence)>26 & strcmp(test_sequence(1:26),'GGGAACGACUCGAGUAGAGUCGAAAA')
        BLANK_OUT5 = 26;
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
function good_idx = figure_out_idx_for_normalization( total_coverage );
COVERAGE_CUTOFF = 10000; % stringent cutoff.
good_idx = find( total_coverage > COVERAGE_CUTOFF );
if length(good_idx) < 10
    COVERAGE_CUTOFF = 1000; % less stringent cutoff.
    good_idx = find( total_coverage > COVERAGE_CUTOFF );
end
if length(good_idx) < 10
    COVERAGE_CUTOFF = 100; % not stringent
    good_idx = find( total_coverage > COVERAGE_CUTOFF );
end
if length(good_idx) < 10
    COVERAGE_CUTOFF = 0; % not stringent
    good_idx = [1:length(total_coverage)];
end
fprintf( 'For normalization, using %d sequences that pass a total coverage cutoff of %d \n',length(good_idx),COVERAGE_CUTOFF)




