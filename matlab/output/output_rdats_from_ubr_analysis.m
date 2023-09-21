function output_rdats_from_ubr_analysis( d, out_prefix, name, output_idx, comments, annotations, condition_comments, condition_annotations, chunk_size, signal_to_noise_cutoff, reads_cutoff, output_bad_profiles);
% output_rdats_from_ubr_analysis( d, out_prefix, name, comments, annotations );
% output_rdats_from_ubr_analysis( d, out_prefix, name, output_idx, comments, annotations, condition_comments, condition_annotations, chunk_size, signal_to_noise_cutoff, reads_cutoff, output_bad_profiles);
%
% Output background-subtracted data from UBR analysis to RDAT format,
%  for profiles that surpass coverage, reads, and signal/noise criteria.
% One RDAT per condition.
% Note that all profiles will be outputted but 'bad' profiles will be
%  marked with the tag "warning:badQuality"
%
% Inputs
%  d      = MATLAB struct() that is output from quick_look_ubr(). 
%            Should contain fields BLANK_OUT5, BLANK_OUT3, r_norm,
%            r_norm_err, sequences, structures. 
%  out_prefix = (string) prefix of RDAT filename, e.g. 'data_sets/release/RH_OpenKnot1'
%  name   = (string) name of library, e.g., 'OpenKnotRound1'.
%  output_idx = indices of which designs to output. Give [] for all. 
%  comments   = comments to put in RDAT
%  annotations = cell of tags -- any annotations to include in ANNOTATION. Make sure to
%                 include tags like {'modifier:1M7','chemical:HEPES:50mM(pH8.0)','chemical:MgCl2:10mM','temperature:24C'}
%
% Optional Inputs
%  condition_comments = cell of Nconditions cells of extra comments for
%                           each condition. 
%  condition_annotations = cell of Nconditions cells with extra annotations
%                           for each conditions
%  chunk_size = max number of sequences per RDAT file (default 10000).
%  signal_to_noise_cutoff = Minimum signal-to-noise (default 1.0).
%  reads_cutoff = Minimum reads (default 100).
%  output_bad_profiles = Output all profiles, including ones that do not
%      pass signal_to_noise and reads cutoffs (which will be marked by
%      'warning:badQuality' in DATA_ANNOTATION). Default 1.
%
% (C) R. Das, Stanford & HHMI, 2023.

Ndesigns    = size(d.r_norm,1);
Nconditions = size(d.r_norm,3);

if ~exist( 'output_idx', 'var') | isempty(output_idx); output_idx = [1:Ndesigns]; end;
if ~exist( 'condition_comments', 'var') | length(condition_comments)==0
    condition_comments = [];
    for i = 1:Nconditions; condition_comments{i} = {}; end;
end
if ~exist( 'condition_annotations', 'var') | length(condition_annotations)==0
    condition_annotations = [];
    for i = 1:Nconditions; condition_annotations{i} = {}; end;
end
if ~exist( 'chunk_size', 'var') | chunk_size == 0; chunk_size = 10000; end;
if ~exist( 'signal_to_noise_cutoff', 'var'); signal_to_noise_cutoff = 1; end
if ~exist( 'reads_cutoff', 'var'); reads_cutoff = 100; end
if length(condition_annotations) ~= Nconditions
    fprintf('Length of condition_annotations %d needs to equal the number of conditions %d.\n',length(condition_annotations),Nconditions)
end
if ~exist('output_bad_profiles','var'); output_bad_profiles = 1; end;
if ~isfield(d,'ids'); d.ids = []; d.titles = []; d.authors = []; end;

if ~exist('RDAT','dir'); mkdir('RDAT'); end;
for i = 1:Nconditions
    good_idx = find( d.signal_to_noise(:,i) >= signal_to_noise_cutoff & d.reads(:,i) >= reads_cutoff);
    fprintf( 'Number of output designs passing cutoff: %3.1f %% (%d out of %d)\n',100*length(intersect(output_idx,good_idx))/length(output_idx),length(intersect(output_idx,good_idx)),length(output_idx));

    annotations_out = [annotations,condition_annotations{i}];
    comments_out = [comments,condition_comments{i}];

    if isfield(d,'norm_val') & length(d.norm_val)>=i; annotations_out = [annotations_out,{sprintf('processing:normalization:value:%.4f',d.norm_val(i))}]; end;
    annotations_out = [annotations_out,{sprintf('processing:UBR-v%s',get_ubr_version())}];

    extra_data_annotations = get_extra_data_annotations_for_eterna( d.ids, d.titles, d.authors, [], d.headers);

    bad_idx = setdiff([1:Ndesigns],good_idx);
    for n = bad_idx; extra_data_annotations{n} = [extra_data_annotations{n},{'warning:badQuality'}]; end

    save_idx = [1:Ndesigns]';
    if ~output_bad_profiles; save_idx = good_idx; end

    filename = sprintf('RDAT/%s_%s.rdat',out_prefix,d.conditions{i});

    num_chunks = floor((length(output_idx)-1)/chunk_size)+1;
    for q = 1:num_chunks        
        chunk_idx = output_idx( ((q-1)*chunk_size+1):min(q*chunk_size,length(output_idx)) );
        save_chunk_idx = intersect(chunk_idx, save_idx);

        filename_chunk = filename;
        if (num_chunks > 1) filename_chunk = strrep(filename,'.rdat',sprintf('_%05d.rdat',q) ); end;

        rdat = output_rdat_from_ubr_analysis( filename_chunk, name, save_chunk_idx, ...
            d.r_norm(:,:,i), d.r_norm_err(:,:,i), ...
            d.sequences, d.structures, d.BLANK_OUT5, d.BLANK_OUT3, ...
            comments_out, annotations_out, extra_data_annotations, d.reads(:,i) );
        %show_rdat(filename_chunk); % check that we can read back in from disk.
    end
end