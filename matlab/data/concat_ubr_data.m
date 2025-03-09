function d_out = concat_ubr_data( d1, d2 )
% d_out = concat_ubr_data( all_d )
% d_out = concat_ubr_data( d1, d2 )
%
% Concatenates UBR dataspaces by sequence.
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
% (C) R. Das, Stanford University, HHMI, 2024

if (nargin == 0), help(mfilename); end
if (nargin == 1), all_d = d1; end
if (nargin == 2), all_d = {d1,d2}; end
assert( iscell(all_d) );
assert( length(all_d)>= 1);


d_out = all_d{1};
d = d_out;
d_out = rmfield(d_out,{'r_norm','r_norm_err','mut_rate_matrix','rfcount_mut_rate_profiles'});

all_idx = {}; 
count = 0;
for q = 1:length( all_d )
    d = all_d{q};        
    if isempty(d); warning( 'Missing data!' ); continue; end;
    all_idx{q} = count+[1:size(d.r_norm,1)];
    count = count + size(d.r_norm,1);
end

% preallocate memory for r and r_err:
r     = zeros(count,size(d.r_norm,2),size(d.r_norm,3));
r_err = zeros(count,size(d.r_norm,2),size(d.r_norm,3));
d_out.reads = zeros(count,size(d.reads,2));
d_out.coverage = zeros(count,size(d.coverage,2));
d_out.signal_to_noise = zeros(count,size(d.signal_to_noise,2));
d_out.sequences = cell(1,count);
d_out.headers   = cell(1,count);
d_out.ids = zeros(1,count)
d_out.titles   = cell(1,count);
d_out.authors   = cell(1,count);
d_out.id_strings   = cell(1,count);

if isfield(d,'coverage_matrix')
    mut_count = zeros('like',d.mut_rate_matrix);
    rfcount   = zeros('like',d.rfcount_mut_rate_profiles);
    coverage_matrix = zeros('like',rfcount);
end

for q = 1:length( all_d )
    fprintf('Concatenating %d of %d\n',q,length(all_d));
    d = all_d{q};
    if isempty(d); warning( 'Missing data!' ); continue; end;
    assert( all(strcmp(d_out.conditions,d.conditions)) );
    assert( all(strcmp(d_out.tags,d.tags)) );
    idx = all_idx{q};
    d_out.reads(idx,:) = d.reads;
    d_out.coverage(idx,:) = d.coverage;
    d_out.signal_to_noise(idx,:) = d.signal_to_noise;
    d_out.sequences(idx) = d.sequences;
    if length(d.headers) == length(idx); d_out.headers(idx) = d.headers; else; warning('Could not create headers since it is not filled correctly.'); end;
    if length(d.ids) == length(idx); d_out.ids(idx) = d.ids; else; warning('Could not create ids since it is not filled correctly.'); end;
    if length(d.authors) == length(idx); d_out.authors(idx) = d.authors; else; warning('Could not create authors since it is not filled correctly.'); end;
    if length(d.titles) == length(idx); d_out.titles(idx) = d.titles; else; warning('Could not create titles since it is not filled correctly.'); end;
    if length(d.id_strings) == length(idx); d_out.id_strings(idx) = d.id_strings; else; warning('Could not create id_strings since it is not filled correctly.'); end;

    r_tmp = []; r_tmp_err = [];
    for i = 1:size(d.r_norm,3)
        r_tmp(:,:,i)     = d.r_norm(:,:,i)*d.norm_val(i);
        r_tmp_err(:,:,i) = d.r_norm_err(:,:,i)*d.norm_val(i);
    end

    r(all_idx{q},:,:) = r_tmp;
    r_err(all_idx{q},:,:) = r_tmp_err;

    if isfield(d,'coverage_matrix');
        mut_count_tmp = [];
        for i =  1:size(d.coverage_matrix,2)
            coverage_matrix_tmp = repmat(reshape( d.coverage_matrix(:,i),[],1), 1, size(d.mut_rate_matrix,2));
            mut_count_tmp(:,:,i) = d.mut_rate_matrix(:,:,i) .* coverage_matrix_tmp;
        end
        mut_count_tmp( isnan(mut_count_tmp) ) = 0;
        mut_count = mut_count + mut_count_tmp;

        coverage_matrix = coverage_matrix + d.coverage_matrix;
        rfcount_tmp = d.rfcount_mut_rate_profiles .* d.coverage_matrix;
        rfcount_tmp( isnan(rfcount_tmp) ) = 0;
        rfcount = rfcount + rfcount_tmp;
    end
    if isfield(d,'rsep_norm');
        for i = 1:size(d.rsep_norm,4)
            rsep(idx,:,:,i)     = d.rsep_norm(:,:,:,i)*d.norm_val(i);
        end
    end
end

d_out.total_coverage = sum(d_out.coverage,2);
norm_idx = figure_out_idx_for_normalization( d_out.total_coverage );
if length(norm_idx)>0; 
    [d_out.r_norm, d_out.r_norm_err,~,d_out.norm_val] = normalize_reactivity(r,r_err,norm_idx,d_out.BLANK_OUT5, d_out.BLANK_OUT3, d_out.conditions, [], d_out.sequences );
end
d_out.norm_idx = norm_idx;

if isfield(d,'coverage_matrix');
    for i =  1:size(d.coverage_matrix,2)
        cov_matrix_repmat = repmat(coverage_matrix(:,i), 1, size(d.mut_rate_matrix,2));
        d_out.mut_rate_matrix(:,:,i) = mut_count(:,:,i)./cov_matrix_repmat;
    end
    d_out.rfcount_mut_rate_profiles = rfcount ./coverage_matrix;
    d_out.coverage_matrix = coverage_matrix;
end

if isfield(d,'rsep_norm');
    d_out.rsep_norm = [];
    for i = 1:size(rsep,4)
        d_out.rsep_norm(:,:,:,i) = rsep(:,:,:,i)/d_out.norm_val(i);
    end
end

