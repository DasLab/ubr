function output_heatmaps_with_mfe_structures_for_design( image_dir, good_idx, r_norm, mfe_structure_map, mfe_structures, mfe_tags, pkg_sort_idx, headers, sequences, ids, titles, authors, BLANK_OUT5, BLANK_OUT3, tags)
% output_heatmaps_with_mfe_structures_for_design( image_dir, good_idx, r_norm, mfe_structure_map, mfe_structures, mfe_tags, pkg_sort_idx, headers, sequences, ids, titles, authors, BLANK_OUT5, BLANK_OUT3, tags)
%
% Deprecated in favor of make_heatmap_with_structures_for_design in OpenKnotScore repo.
%
% (C) R. Das, Stanford/HHMI 2023.

if ~exist('tags','var') tags = {'SHAPE, no Mg2+','SHAPE, +Mg2+'}; end
if ~exist('create_figure_window','var') create_figure_window = 1; end

if ~exist( 'output_heatmaps_with_structures_for_design', 'file' )
  fprintf( 'ERROR! You need to install OpenKnotScore from https://github.com/eternagame/OpenKnotScore, and add OpenKnotScore/matlab and subdirectories to your path.' );
  return;
end

output_heatmaps_with_structures_for_design( image_dir, good_idx, r_norm, mfe_structure_map, mfe_structures, mfe_tags, pkg_sort_idx, headers, sequences, ids, titles, authors, BLANK_OUT5, BLANK_OUT3, tags)

fprintf( 'WARNING! output_heatmaps_with_mfe_structures_for_design has been deprecated in favor of output_heatmaps_with_structures_for_design in the OpenKnotScore repository at https://github.com/eternagame/OpenKnotScore\n' );
