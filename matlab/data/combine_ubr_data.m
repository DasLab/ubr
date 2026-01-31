function d_out = combine_ubr_data( d1, d2, RENORM )
% d_out = combine_ubr_data( all_d )
% d_out = combine_ubr_data( d1, d2 )
%
% Combines two data sets for the same sequences but in different conditions
%  Also works on data sets with 2d r, r_norm data.
%
%  See also:
%     MERGE_UBR_DATA to actually take averages.
%     CONCAT_UBR_DATA to concatenate down sequence dimension rather than
%     across conditions
%
% Inputs
%  all_d = cell of MATLAB data structures with fields like r_norm, reads, etc. (as
%        output from quick_look_ubr),
%              OR
%  d1 = MATLAB data structure 1
%  d2 = MATLAB data structure 2
%
% Output
%  d_out = MATLAB data structure with merged
%
% (C) R. Das, Stanford University, HHMI, 2025

if (nargin == 0), help(mfilename); end
if (nargin == 1), all_d = d1; end
if (nargin == 2 & isstruct(d2) ), all_d = {d1,d2}; end

assert( iscell(all_d) );
assert( length(all_d)>= 1);

d_out = all_d{1};

for q = 2:length( all_d )
    num_cidx = length(d_out.conditions);
    assert( length(all_d{q}.sequences) == length(d_out.sequences) );
    assert( all(strcmp(all_d{q}.sequences,d_out.sequences)) );
    NDIMS = length(size(d_out.r_norm));
    if size(d_out.r_norm,NDIMS)~=num_cidx & size(d_out.r_norm,NDIMS+1)==num_cidx; NDIMS = NDIMS+1; end;
    d_out.r_norm = cat(NDIMS, d_out.r_norm, all_d{q}.r_norm);
    d_out.r_norm_err = cat(NDIMS, d_out.r_norm_err, all_d{q}.r_norm_err);
    d_out.conditions =[d_out.conditions, all_d{q}.conditions];
    d_out.reads = cat(2, d_out.reads, all_d{q}.reads);
    d_out.signal_to_noise = cat(2, d_out.signal_to_noise, all_d{q}.signal_to_noise);

    if isfield(d_out,'coverage'); 
        d_out.coverage = cat(2, d_out.coverage, all_d{q}.coverage);
        d_out.total_coverage = d_out.total_coverage + all_d{q}.total_coverage;
    end;

    if ~isfield(d_out,'tags') continue; end
    num_tags = length(d_out.tags);
    num_shape_nomod_idx = length(d_out.shape_nomod_idx);
    shape_nomod_idx = all_d{q}.shape_nomod_idx;
    for k = 1:length(shape_nomod_idx); 
        d_out.shape_nomod_idx{ num_shape_nomod_idx + k } = shape_nomod_idx{k} + num_tags; 
    end;
    d_out.tags = [d_out.tags, all_d{q}.tags];

    if ~isfield(d_out,'norm_val') continue; end
    d_out.norm_val = [d_out.norm_val, all_d{q}.norm_val];    
    d_out.norm_idx = [];  % will fix in a minute!

    if ~isfield(d_out,'mut_rate_matrix') continue; end
    d_out.mut_rate_matrix = cat(3, d_out.mut_rate_matrix, all_d{q}.mut_rate_matrix);
    d_out.rfcount_mut_rate_profiles = cat(2, d_out.rfcount_mut_rate_profiles, all_d{q}.rfcount_mut_rate_profiles);
    d_out.coverage_matrix = cat(2, d_out.coverage_matrix, all_d{q}.coverage_matrix);
end


if ~isfield(d_out,'norm_val') return; end

% RENORM
norm_idx = figure_out_idx_for_normalization( d_out.total_coverage );
if ~isempty(norm_idx)
    for i = 1:size(d_out.r_norm,3)
        r(:,:,i)     = d_out.r_norm(:,:,i)*d_out.norm_val(i);
        r_err(:,:,i) = d_out.r_norm_err(:,:,i)*d_out.norm_val(i);
    end
    [d_out.r_norm, d_out.r_norm_err,~,d_out.norm_val] = normalize_reactivity(r,r_err,norm_idx,d_out.BLANK_OUT5, d_out.BLANK_OUT3, d_out.conditions, [], d_out.sequences );
end

d_out.norm_idx = norm_idx;


