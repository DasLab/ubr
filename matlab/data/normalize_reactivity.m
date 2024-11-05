function [r, r_err,r_nomod, norm_val] = normalize_reactivity(r,r_err,good_idx,BLANK_OUT5, BLANK_OUT3, tags_conditions, r_nomod, sequences );
% [r, r_err,r_nomod, norm_val]  = normalize_reactivity(r,r_err,good_idx,BLANK_OUT5, BLANK_OUT3, tags_conditions);
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
%  r       = [Ndesign x Nres] Reactivity matrix, normalized.
%  r_err   = [Ndesign x Nres] Reactivity error matrix, normalized.
%  r_nomod = [Ndesign x Nres] Nomod value matrix, normalized, can be
%                      used for systematic error estimate. 
%                       (Zeros if r_nomod not provided)
%  norm_val = values (in units of fraction modified) used for normalization
%                      
%
% (C) R. Das, HHMI/Stanford University 2023.

r = single(r);
r_err = single(r_err);

norm_val = [];
if ~exist( 'tags_conditions','var'); tags_conditions = repmat({''},1,size(r,3)); end;
if ~exist('r_nomod','var'); r_nomod = []; end

% "Blank out" 5' and 3' sequence as NaN.
r(:,1:BLANK_OUT5,:) = NaN;
r(:,(end-BLANK_OUT3+1):end,:) = NaN;
if exist( 'sequences','var') & length(sequences)>0
    unique_lens = unique(cellfun(@length,sequences));
    if length(unique_lens > 1)
        for i = 1:length(sequences);
            N = length(sequences{i});
            r(i,(N-BLANK_OUT3+1):end,:) = NaN;
        end
    else
        N = unique_lens;
        r(:,(N-BLANK_OUT3+1):end,:) = NaN;
    end
end

% Get 90th percentile values over the sequences with norm_idx
for i = 1:size(r,3)
    vals = r(good_idx,:,i);
    vals = vals( find(~isnan(vals)));
    vals_sort = sort(vals(:));
    val_norm = 0.0;
    if length(vals_sort)>0; 
        val_norm = vals_sort( round(0.90 * length(vals_sort) )); 
    end;
    % If 90th percentile of all vals is 0.0, look only at non-zero vals
    % (and issue warning!)
    if val_norm==0.0;     
        vals = vals(vals~=0.0);  
        vals_sort = sort(vals(:)); 
        if length(vals_sort)>0; 
            val_norm = vals_sort( round(0.90 * length(vals_sort) )); 
            warning(sprintf('Had to expand number of profiles to all for normalization of %d %s',i,tags_conditions{i}));  
        end;
    end
    fprintf( 'Normalizing reactivity profiles for %d %s with value %f.\n', i, tags_conditions{i}, val_norm);
    norm_val(i) = val_norm;
end;

% Need to normalize 'in place', and in chunks to avoid memory issues
chunk_size = 1e5;
Nseq = size(r,1);
num_chunks = ceil(Nseq/chunk_size);
for q = 1:num_chunks
    chunk_start = 1+(q-1)*chunk_size;
    chunk_end = min(chunk_start + chunk_size-1,Nseq);
    seq_range = [chunk_start:chunk_end];

    % use tensor product for speed -- looks like an issue though where
    % tensorprod seems to make things nan?
    % r(seq_range,:,:)     = tensorprod(r(seq_range,:,:),    diag(1./norm_val),3,2);
    % r_err(seq_range,:,:) = tensorprod(r_err(seq_range,:,:),diag(1./norm_val),3,2);
    norm_val_rep = repmat(permute(norm_val,[3 1 2]),length(seq_range),size(r,2));   
    r(seq_range,:,:)     = r(seq_range,:,:)./norm_val_rep;
    r_err(seq_range,:,:) = r_err(seq_range,:,:)./norm_val_rep;
    if exist('r_nomod','var') & ~isempty(r_nomod); 
        %r_nomod(seq_range,:,:) = tensorprod(r_nomod(seq_range,:,:),diag(1./norm_val),3,1);
        r_nomod(seq_range,:,:) = r_nomod(seq_range,:,:)./norm_val_rep;
    end;
end
