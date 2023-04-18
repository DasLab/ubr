function [r_norm, r_norm_err] = normalize_reactivity(r,r_err,sequences,good_idx,BLANK_OUT5, BLANK_OUT3, tags_conditions );
% [r_norm, r_norm_err] = normalize_reactivity(r,r_err,good_idx,BLANK_OUT5, BLANK_OUT3, tags_conditions);
%
% Inputs
%  r     = [Ndesign x Nres x Nconditions] Reactivity matrix, as fraction at each position that leads to mutation.
%  r_err = [Ndesign x Nres x Nconditions] Error of reactivity matrix, as fraction at each position that leads to mutation.
%  good_idx = [list of integers] index of designs to use for normalization
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues
%  tags_conditions = tags for each condition in r (leave out to show
%                           blanks)
%
% Outputs
%  r_norm     = [Ndesign x Nres] Reactivity matrix, normalized.
%  r_norm_err = [Ndesign x Nres] Reactivity error matrix, normalized.
%
% (C) R. Das, HHMI/Stanford University 2023.

N = size(r,2);
which_pos = [(BLANK_OUT5+1):(N-BLANK_OUT3)];
if ~exist( 'tags_conditions','var'); tags_conditions = repmat({''},1,size(r,3)); end;
nt_types = {'A','C','G','U'};

for i = 1:size(r,3)
    fprintf( '\n');        
    for nt_type_idx = 1:4;
        nt_type = nt_types{nt_type_idx};
        vals = [];
        for k = good_idx';
            good_pos = intersect(which_pos, strfind(sequences{k},nt_type));
            vals = [vals, r(k,good_pos,i)];
        end
        vals_sort = sort(vals(:));
        val_norm = vals_sort( round(0.90 * length(vals_sort) ));
        fprintf( 'Normalizing reactivity profiles for %d %s at %s with value %f.\n', i, tags_conditions{i}, nt_type, val_norm);
        for k = 1:size(r,1)
            good_pos = strfind(sequences{k},nt_type);
            r_norm(k,good_pos,i) = r(k,good_pos,i)/val_norm;
            r_norm_err(k,good_pos,i) = r_err(k,good_pos,i)/val_norm;
        end
    end
end
