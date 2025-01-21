function d_out = concat_ubr_data( d, slice_idx )
% d_out = concat_ubr_data( d, slice_idx )
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


