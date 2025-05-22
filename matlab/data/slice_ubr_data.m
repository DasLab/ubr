function d_out = slice_ubr_data( d, slice_idx )
% d_out = slice_ubr_data( d, slice_idx )
%
% Slices UBR dataspaces to just sequences specified by slice_idx
%
% Inputs
%  d = MATLAB data structure
%  slice_idx = vector of integers, which sequences to slice out
%
% Output
%  d_out = MATLAB data structure with merged
%
% (C) R. Das, Stanford University, HHMI, 2025

if (nargin < 2), help(mfilename); end

d_out.sequences = d.sequences(slice_idx);
d_out.conditions = d.conditions;
d_out.reads = d.reads( slice_idx,: );
d_out.signal_to_noise = d.signal_to_noise( slice_idx, :);
d_out.BLANK_OUT5 = d.BLANK_OUT5;
d_out.BLANK_OUT3 = d.BLANK_OUT3;
d_out.headers = d.headers( slice_idx );
% d_out.authors = d.authors( slice_idx );
% d_out.ids = d.ids( slice_idx );
% d_out.titles = d.titles( slice_idx );
% d_out.id_strings = d.id_strings( slice_idx );
d_out.shape_nomod_idx = d.shape_nomod_idx;
d_out.tags = d.tags;
d_out.coverage = d.coverage( slice_idx, :);
d_out.total_coverage = d.total_coverage( slice_idx, :);
d_out.norm_idx = []; % necessary incomplete
d_out.norm_val = d.norm_val;
d_out.structures = {};
if ~isempty( d.structures ) & length(d.structures) > max(slice_idx);  d_out.structures = d.structure(slice_idx); end;
d_out.structure_map = [];
if ~isempty( d.structure_map ) & size(d.structure_map,1) > max(slice_idx);  d_out.structure_map = d.structure_map(slice_idx,:); end;
d_out.filedir = d.filedir;
d_out.structure_csv_file = d.structure_csv_file;
d_out.coverage_matrix = d.coverage_matrix;
d_out.r_norm = d.r_norm( slice_idx, :, :);
d_out.r_norm_err = d.r_norm_err( slice_idx, :, :);
d_out.mut_rate_matrix = d.mut_rate_matrix;
d_out.rfcount_mut_rate_profiles = d.rfcount_mut_rate_profiles;
if isfield(d,'rsep_norm'); d_out.rsep_norm = d.rsep_norm( slice_idx, :, :, :); end;

r = []; r_err = [];
for i = 1:size(d.r_norm,3)
    r(:,:,i)     = d_out.r_norm(:,:,i)*d.norm_val(i);
    r_err(:,:,i) = d_out.r_norm_err(:,:,i)*d.norm_val(i);
    if isfield(d,'rsep_norm'); rsep(:,:,:,i) = d_out.rsep_norm(:,:,:,i)*d.norm_val(i); end;
end

d_out.total_coverage = sum(d_out.coverage,2);
norm_idx = figure_out_idx_for_normalization( d_out.total_coverage );
if length(norm_idx)>0; 
    [d_out.r_norm, d_out.r_norm_err,~,d_out.norm_val] = normalize_reactivity(r,r_err,norm_idx,d_out.BLANK_OUT5, d_out.BLANK_OUT3, d_out.conditions, [], d_out.sequences );
end
d_out.norm_idx = norm_idx;

if isfield(d,'rsep_norm');
    d_out.rsep_norm = [];
    for i = 1:size(rsep,4)
        d_out.rsep_norm(:,:,:,i) = rsep(:,:,:,i)/d_out.norm_val(i);
    end
end
