function [r_norm, r_norm_err,r_norm_nomod,norm_val] = normalize_reactivity(r,r_err,good_idx,BLANK_OUT5, BLANK_OUT3, tags_conditions, r_nomod, sequences );
% [r_norm, r_norm_err] = normalize_reactivity(r,r_err,good_idx,BLANK_OUT5, BLANK_OUT3, tags_conditions);
%
% Inputs
%  r     = [Ndesign x Nres x Nconditions] Reactivity matrix, as fraction at each position that leads to mutation.
%  r_err = [Ndesign x Nres x Nconditions] Error of reactivity matrix, as fraction at each position that leads to mutation.
%  good_idx = [list of integers] index of designs to use for normalization
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%  tags_conditions = [cell of Nconditions strings] tags for each condition in r (leave out to show
%                           blanks)
%  r_nomod = [Optional] [Ndesign x Nres x Nmodconditions] No mod values subtracted out to get r
%  sequences = [Optional] cell of Ndesign sequences, used to properly blank
%                  out 5' and 3' when sequences are variable length.
%
% Outputs
%  r_norm     = [Ndesign x Nres] Reactivity matrix, normalized.
%  r_norm_err = [Ndesign x Nres] Reactivity error matrix, normalized.
%  r_norm_nomod = [Ndesign x Nres] Nomod value matrix, normalized, can be
%                      used for systematic error estimate. 
%                       (Zeros if r_nomod not provided)
%  norm_val = values (in units of fraction modified) used for normalization
%                      
%
% (C) R. Das, HHMI/Stanford University 2023.

r_norm = [];
r_norm_err = [];
r_norm_nomod = [];
norm_val = [];
if ~exist( 'tags_conditions','var'); tags_conditions = repmat({''},1,size(r,3)); end;

% "Blank out" 5' and 3' sequence as NaN.
r(:,1:BLANK_OUT5,:) = NaN;
r(:,(end-BLANK_OUT3+1):end,:) = NaN;
if exist( 'sequences','var') & length(sequences)>0
    for i = 1:length(sequences);
        N = length(sequences{i});
        r(:,(N-BLANK_OUT3+1):end,:) = NaN;
    end
end
for i = 1:size(r,3)
    vals = r(good_idx,:,i);
    vals = vals( find(~isnan(vals)));
    vals_sort = sort(vals(:));
    val_norm = vals_sort( round(0.90 * length(vals_sort) ));
    fprintf( 'Normalizing reactivity profiles for %d %s with value %f.\n', i, tags_conditions{i}, val_norm);
    r_norm(:,:,i) = single(r(:,:,i))/single(val_norm);
    r_norm_err(:,:,i) = single(r_err(:,:,i))/single(val_norm);
    r_norm_nomod(:,:,i) = 0*single(r_norm(:,:,i));
    if exist('r_nomod','var') & ~isempty(r_nomod); r_norm_nomod(:,:,i) = single(r_nomod(:,:,i))/single(val_norm); end;
    norm_val(i) = val_norm;
end;
